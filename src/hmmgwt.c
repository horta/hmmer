#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gencode.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

#define RAISE_EWRITE ESL_EXCEPTION_SYS(eslEWRITE, "write failed")
#define XRAISE_EWRITE ESL_XEXCEPTION_SYS(eslEWRITE, "write failed")

#define FPRINTF(stream, format, ...)                                                   \
  do {                                                                                 \
    if (fprintf(stream, format, __VA_ARGS__) < 0)                                      \
      ESL_EXCEPTION_SYS(eslEWRITE, "write failed");                                    \
  } while (0)

#define PUTS(str)                                                                      \
  do {                                                                                 \
    if (puts(str) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");                   \
  } while (0)

/* struct cfg_s : "Global" application configuration shared by all
 * threads/processes
 *
 * This structure is passed to routines within main.c, as a means of
 * semi-encapsulation of shared data amongst different parallel processes
 * (threads).
 */
struct cfg_s {
  const ESL_GETOPTS *go;
  char *seq_fpath;
  char *hmm_fpath;
};

ESL_GETOPTS *create_getopts(int argc, char **argv);
int process_cmdline(const ESL_GETOPTS *go, char **ret_hmm_filepath,
                    char **ret_seq_filepath);
int serial_master(const ESL_GETOPTS *go, struct cfg_s *cfg);

int main(int argc, char **argv)
{
  impl_Init();
  p7_FLogsumInit();

  ESL_GETOPTS *go = create_getopts(argc, argv);

  struct cfg_s cfg = {.go = go, .seq_fpath = NULL, .hmm_fpath = NULL};

  int status = process_cmdline(go, &cfg.hmm_fpath, &cfg.seq_fpath);
  if (status) return status;

  status = serial_master(go, &cfg);
  if (status) return status;

  esl_getopts_Destroy(go);
  return status;
}

void usage_err(const ESL_GETOPTS *go);
void show_help(const ESL_GETOPTS *go);

int process_cmdline(const ESL_GETOPTS *go, char **ret_hmm_filepath,
                    char **ret_seq_filepath)
{
  if (esl_opt_GetBoolean(go, "-h") == TRUE) show_help(go);

  if (esl_opt_ArgNumber(go) != 2) {
    PUTS("Incorrect number of command line arguments.");
    usage_err(go);
  }

  if ((*ret_hmm_filepath = esl_opt_GetArg(go, 1)) == NULL) {
    PUTS("Failed to get <hmmfile> argument on command line.");
    usage_err(go);
  }

  if ((*ret_seq_filepath = esl_opt_GetArg(go, 2)) == NULL) {
    PUTS("Failed to get <seqfile> argument on command line.");
    usage_err(go);
  }

  return eslOK;
}

struct GenTable {
  ESL_ALPHABET *dna;
  ESL_ALPHABET *amino;
  ESL_GENCODE *gcode;
};

struct GenTable create_gentable();
void destroy_gentable(struct GenTable gentbl);
ESL_SQFILE *open_seqfile(const char *filepath, int fmt);
int get_seqfile_format(const ESL_GETOPTS *go);
int loop_query_seq(const struct cfg_s *cfg, ESL_SQFILE *seqfile,
                   ESL_SQ *sq, struct GenTable gentbl);
void load_hmm(const char *fpath, ESL_ALPHABET **ret_abc, P7_HMM **ret_hmm);
static int output_header(FILE *ofp, const ESL_GETOPTS *go, const char *hmm_fpath,
                         const char *seq_fpath);

int serial_master(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int status = eslOK;
  output_header(stdout, go, cfg->hmm_fpath, cfg->seq_fpath);

  struct GenTable gentbl = create_gentable();

  ESL_ALPHABET *abc = NULL;
  P7_HMM *hmm = NULL;
  load_hmm(cfg->hmm_fpath, &abc, &hmm);

  ESL_SQFILE *seqfile = open_seqfile(cfg->seq_fpath, get_seqfile_format(go));
  esl_sqfile_SetDigital(seqfile, abc);
  ESL_SQ *sq = esl_sq_CreateDigital(abc);

  loop_query_seq(cfg, seqfile, sq, gentbl);

  esl_sqfile_Close(seqfile);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  destroy_gentable(gentbl);
  return status;
}

struct GenTable create_gentable()
{
  struct GenTable gentbl;
  gentbl.dna = esl_alphabet_Create(eslDNA);
  if (!gentbl.dna) esl_fatal("Could not create the DNA alphabet");

  gentbl.amino = esl_alphabet_Create(eslAMINO);
  if (!gentbl.amino) esl_fatal("Could not create the AMINO alphabet");

  gentbl.gcode = esl_gencode_Create(gentbl.dna, gentbl.amino);
  if (!gentbl.gcode) esl_fatal("Could not create the genetic code");

  return gentbl;
}

void destroy_gentable(struct GenTable gentbl)
{
  esl_alphabet_Destroy(gentbl.dna);
  esl_alphabet_Destroy(gentbl.amino);
  esl_gencode_Destroy(gentbl.gcode);
}

struct WORKER_INFO {
  P7_BG *bg;
  P7_PIPELINE *pli;
  P7_TOPHITS *th;
  P7_OPROFILE *om;
  ESL_SQ *sq;
};

int serial_loop(struct WORKER_INFO *info, P7_HMMFILE *hfp);

void seq_fail(const char *filepath, int status, ESL_SQFILE *seqfile);

void workinfo_setup(struct WORKER_INFO *info, const ESL_GETOPTS *go, P7_HMMFILE *hmmfile,
                    ESL_SQ *sq)
{
  info->th = p7_tophits_Create();
  info->pli = p7_pipeline_Create(go, 100, 100, FALSE, p7_SCAN_MODELS);
  info->pli->hfp = hmmfile;

  p7_pli_NewSeq(info->pli, sq);
  info->sq = sq;
}

int loop_query_seq(const struct cfg_s *cfg, ESL_SQFILE *seqfile,
                   ESL_SQ *sq, struct GenTable gentbl)
{
  P7_HMMFILE *hmmfile = NULL;

  int hstatus = eslOK;
  int sstatus = eslOK;
  while ((sstatus = esl_sqio_Read(seqfile, sq)) == eslOK) {
    hstatus = p7_hmmfile_OpenE(cfg->hmm_fpath, p7_HMMDBENV, &hmmfile, NULL);
    if (hstatus)
      p7_Fail("Unexpected error %d in opening hmm file %s.\n", hstatus, cfg->hmm_fpath);

    FPRINTF(stdout, "Query:       %s  [L=%ld]\n", sq->name, (long)sq->n);
    if (sq->acc[0]) FPRINTF(stdout, "Accession:   %s\n", sq->acc);
    if (sq->desc[0]) FPRINTF(stdout, "Description: %s\n", sq->desc);

    struct WORKER_INFO info;
    workinfo_setup(&info, cfg->go, hmmfile, sq);

    p7_hmmfile_Close(hmmfile);
  }
  if (sstatus != eslEOF) seq_fail(cfg->seq_fpath, sstatus, seqfile);

  return eslOK;
}

/* int serial_loop(struct WORKER_INFO *info, P7_HMMFILE *hfp) */
/* { */
/*   int status; */

/*   P7_OPROFILE *om; */
/*   ESL_ALPHABET *abc = NULL; */
/*   /1* Main loop: *1/ */
/*   while ((status = p7_oprofile_ReadMSV(hfp, &abc, &om)) == eslOK) { */
/*     p7_pli_NewModel(info->pli, om, info->bg); */
/*     p7_bg_SetLength(info->bg, info->qsq->n); */
/*     p7_oprofile_ReconfigLength(om, info->qsq->n); */

/*     status = p7_Pipeline(info->pli, om, info->bg, info->qsq, NULL, info->th, NULL);
 */
/*     if (status == eslEINVAL) p7_Fail(info->pli->errbuf); */

/*     p7_oprofile_Destroy(om); */
/*     p7_pipeline_Reuse(info->pli); */
/*   } */

/*   esl_alphabet_Destroy(abc); */

/*   return status; */
/* } */

int get_seqfile_format(const ESL_GETOPTS *go)
{
  int seqfile_fmt = eslSQFILE_UNKNOWN;

  if (esl_opt_IsOn(go, "--tformat")) {
    seqfile_fmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--tformat"));
    if (seqfile_fmt == eslSQFILE_UNKNOWN)
      p7_Fail("%s is not a recognized sequence database file format\n",
              esl_opt_GetString(go, "--tformat"));
  }
  return seqfile_fmt;
}

ESL_SQFILE *open_seqfile(const char *filepath, int fmt)
{
  ESL_SQFILE *seqfile = NULL;

  int status = esl_sqfile_Open(filepath, fmt, p7_SEQDBENV, &seqfile);
  if (status == eslENOTFOUND)
    p7_Fail("Failed to open sequence file %s for reading\n", filepath);
  else if (status == eslEFORMAT)
    p7_Fail("Sequence file %s is empty or misformatted\n", filepath);
  else if (status == eslEINVAL)
    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)
    p7_Fail("Unexpected error %d opening sequence file %s\n", status, filepath);

  return seqfile;
}

void seq_fail(const char *filepath, int status, ESL_SQFILE *seqfile)
{
  if (status == eslEFORMAT)
    esl_fatal("Parse failed (sequence file %s):\n%s\n", filepath,
              esl_sqfile_GetErrorBuf(seqfile));
  else
    esl_fatal("Unexpected error %d reading sequence file %s", status, filepath);
}

void hmm_fail(const char *filepath, int status)
{
  switch (status) {
    case eslEOD:
      p7_Fail("read failed, HMM file %s may be truncated?", filepath);
      break;
    case eslEFORMAT:
      p7_Fail("bad file format in HMM file %s", filepath);
      break;
    case eslEINCOMPAT:
      p7_Fail("HMM file %s contains different alphabets", filepath);
      break;
    case eslEOF: /* do nothing. EOF is what we want. */
      break;
    default:
      p7_Fail("Unexpected error (%d) in reading HMMs from %s", status, filepath);
  }
}

P7_HMMFILE *open_hmmfile(const char *filepath)
{
  P7_HMMFILE *hmmfile = NULL;
  char errbuf[eslERRBUFSIZE];

  int status = p7_hmmfile_OpenE(filepath, NULL, &hmmfile, errbuf);
  if (status == eslENOTFOUND)
    p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n",
            filepath, errbuf);
  else if (status == eslEFORMAT)
    p7_Fail("File format problem in trying to open HMM file %s.\n%s\n", filepath,
            errbuf);
  else if (status != eslOK)
    p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n", status, filepath,
            errbuf);

  return hmmfile;
}

void load_hmm(const char *fpath, ESL_ALPHABET **ret_abc, P7_HMM **ret_hmm)
{
  P7_HMMFILE *hmmfile = open_hmmfile(fpath);
  int hstatus = p7_hmmfile_Read(hmmfile, ret_abc, ret_hmm);
  if (hstatus) hmm_fail(fpath, hstatus);
  p7_hmmfile_Close(hmmfile);
}

static ESL_OPTIONS options[] = {
    {"-h", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL,
     "show brief help on version and usage", 1},
    {"--tformat", eslARG_STRING, NULL, NULL, NULL, NULL, NULL, NULL,
     "assert target <seqfile> is in format <s>: no autodetection", 12},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
};

char usage[] = "[options] <hmmfile> <seqfile>";
char banner[] = "match DNA sequence against protein profile";

int output_header(FILE *ofp, const ESL_GETOPTS *go, const char *hmm_fpath,
                  const char *seq_fpath)
{
  p7_banner(ofp, go->argv[0], banner);
  const char *s = NULL;

  if (fprintf(ofp, "# Protein HMM file:         %s\n", hmm_fpath) < 0) RAISE_EWRITE;
  if (fprintf(ofp, "# DNA sequence file:        %s\n", seq_fpath) < 0) RAISE_EWRITE;

  if (esl_opt_IsUsed(go, "--tformat")) {
    s = esl_opt_GetString(go, "--tformat");
    if (fprintf(ofp, "# targ <seqfile> format asserted:  %s\n", s) < 0) RAISE_EWRITE;
  }

  s = "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n";
  if (fputs(s, ofp) < 0) RAISE_EWRITE;

  return eslOK;
}

void usage_err(const ESL_GETOPTS *go)
{
  int status = eslFAIL;

  esl_usage(stdout, go->argv[0], usage);
  if (puts("\nwhere most common options are:") < 0) XRAISE_EWRITE;

  esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
  if (printf("\nTo see more help on available options, do %s -h.\n\n", go->argv[0]) < 0)
    XRAISE_EWRITE;

  esl_getopts_Destroy(go);
  exit(status);

ERROR:
  esl_getopts_Destroy(go);
  exit(status);
}

void show_help(ESL_GETOPTS *go)
{
  int status = eslOK;
  p7_banner(stdout, go->argv[0], banner);
  esl_usage(stdout, go->argv[0], usage);

  if (puts("\nBasic options:") < 0) XRAISE_EWRITE;
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80);

  if (puts("\nOptions directing output:") < 0) XRAISE_EWRITE;
  esl_opt_DisplayHelp(stdout, go, 2, 2, 80);

  if (puts("\nOptions controlling reporting thresholds:") < 0) XRAISE_EWRITE;
  esl_opt_DisplayHelp(stdout, go, 4, 2, 80);

  if (puts("\nOptions controlling inclusion (significance) thresholds:") < 0)
    XRAISE_EWRITE;
  esl_opt_DisplayHelp(stdout, go, 5, 2, 80);

  if (puts("\nOptions controlling model-specific thresholding:") < 0) XRAISE_EWRITE;
  esl_opt_DisplayHelp(stdout, go, 6, 2, 80);

  if (puts("\nOptions controlling acceleration heuristics:") < 0) XRAISE_EWRITE;
  esl_opt_DisplayHelp(stdout, go, 7, 2, 80);

  if (puts("\nOther expert options:") < 0) XRAISE_EWRITE;
  esl_opt_DisplayHelp(stdout, go, 12, 2, 80);

  if (puts("\nTranslation options:") < 0) XRAISE_EWRITE;
  esl_opt_DisplayHelp(stdout, go, 15, 2, 80);

  exit(eslOK);

ERROR:
  esl_getopts_Destroy(go);
  exit(status);
}

ESL_GETOPTS *create_getopts(int argc, char **argv)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  int status = eslOK;

  if (esl_opt_ProcessEnvironment(go)) {
    if (printf("Failed to process environment: %s\n", go->errbuf) < 0) XRAISE_EWRITE;
    usage_err(go);
  }

  if (esl_opt_ProcessCmdline(go, argc, argv)) {
    if (printf("Failed to parse command line: %s\n", go->errbuf) < 0) XRAISE_EWRITE;
    usage_err(go);
  }

  if (esl_opt_VerifyConfig(go)) {
    if (printf("Failed to parse command line: %s\n", go->errbuf) < 0) XRAISE_EWRITE;
    usage_err(go);
  }

  return go;

ERROR:
  esl_getopts_Destroy(go);
  exit(status);
}
