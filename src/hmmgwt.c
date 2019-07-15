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
#include "hmmgwt.h"

#define REPOPTS "-E,-T,--cut_ga,--cut_nc,--cut_tc"
#define DOMREPOPTS "--domE,--domT,--cut_ga,--cut_nc,--cut_tc"
#define INCOPTS "--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
#define INCDOMOPTS "--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
#define THRESHOPTS                                                                     \
  "-E,-T,--domE,--domT,--incE,--incT,--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"

/* struct cfg_s : "Global" application configuration shared by all
 * threads/processes
 *
 * This structure is passed to routines within main.c, as a means of
 * semi-encapsulation of shared data amongst different parallel processes
 * (threads).
 */
struct cfg_s {
  const ESL_GETOPTS *go;
  char *sfpath;
  char *hfpath;
  int textw;
};

struct WORKER_INFO {
  const struct cfg_s *cfg;
  P7_BG *bg;
  P7_PIPELINE *pli;
  P7_TOPHITS *th;
  P7_OPROFILE *om;
  ESL_SQ *sq;
};

static ESL_OPTIONS options[] = {
    /* name           type          default  env  range toggles  reqs   incomp help
       docgroup*/
    {"-h", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL,
     "show brief help on version and usage", 1},
    /* Control of output */
    {"-o", eslARG_OUTFILE, NULL, NULL, NULL, NULL, NULL, NULL,
     "direct output to file <f>, not stdout", 2},
    {"--tblout", eslARG_OUTFILE, NULL, NULL, NULL, NULL, NULL, NULL,
     "save parseable table of per-sequence hits to file <f>", 2},
    {"--domtblout", eslARG_OUTFILE, NULL, NULL, NULL, NULL, NULL, NULL,
     "save parseable table of per-domain hits to file <f>", 2},
    {"--pfamtblout", eslARG_OUTFILE, NULL, NULL, NULL, NULL, NULL, NULL,
     "save table of hits and domains to file, in Pfam format <f>", 2},
    {"--acc", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL,
     "prefer accessions over names in output", 2},
    {"--noali", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL,
     "don't output alignments, so output is smaller", 2},
    {"--notextw", eslARG_NONE, NULL, NULL, NULL, NULL, NULL, "--textw",
     "unlimit ASCII text output line width", 2},
    {"--textw", eslARG_INT, "120", NULL, "n>=120", NULL, NULL, "--notextw",
     "set max width of ASCII text output lines", 2},
    /* Control of reporting thresholds */
    {"-E", eslARG_REAL, "10.0", NULL, "x>0", NULL, NULL, REPOPTS,
     "report models <= this E-value threshold in output", 4},
    {"-T", eslARG_REAL, FALSE, NULL, NULL, NULL, NULL, REPOPTS,
     "report models >= this score threshold in output", 4},
    {"--domE", eslARG_REAL, "10.0", NULL, "x>0", NULL, NULL, DOMREPOPTS,
     "report domains <= this E-value threshold in output", 4},
    {"--domT", eslARG_REAL, FALSE, NULL, NULL, NULL, NULL, DOMREPOPTS,
     "report domains >= this score cutoff in output", 4},
    /* Control of inclusion (significance) thresholds: */
    {"--incE", eslARG_REAL, "0.01", NULL, "x>0", NULL, NULL, INCOPTS,
     "consider models <= this E-value threshold as significant", 5},
    {"--incT", eslARG_REAL, FALSE, NULL, NULL, NULL, NULL, INCOPTS,
     "consider models >= this score threshold as significant", 5},
    {"--incdomE", eslARG_REAL, "0.01", NULL, "x>0", NULL, NULL, INCDOMOPTS,
     "consider domains <= this E-value threshold as significant", 5},
    {"--incdomT", eslARG_REAL, FALSE, NULL, NULL, NULL, NULL, INCDOMOPTS,
     "consider domains >= this score threshold as significant", 5},
    /* Model-specific thresholding for both reporting and inclusion */
    {"--cut_ga", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, THRESHOPTS,
     "use profile's GA gathering cutoffs to set all thresholding", 6},
    {"--cut_nc", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, THRESHOPTS,
     "use profile's NC noise cutoffs to set all thresholding", 6},
    {"--cut_tc", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, THRESHOPTS,
     "use profile's TC trusted cutoffs to set all thresholding", 6},
    /* Control of acceleration pipeline */
    {"--max", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, "--F1,--F2,--F3",
     "Turn all heuristic filters off (less speed, more power)", 7},
    {"--F1", eslARG_REAL, "0.02", NULL, NULL, NULL, NULL, "--max",
     "MSV threshold: promote hits w/ P <= F1", 7},
    {"--F2", eslARG_REAL, "1e-3", NULL, NULL, NULL, NULL, "--max",
     "Vit threshold: promote hits w/ P <= F2", 7},
    {"--F3", eslARG_REAL, "1e-5", NULL, NULL, NULL, NULL, "--max",
     "Fwd threshold: promote hits w/ P <= F3", 7},
    {"--nobias", eslARG_NONE, NULL, NULL, NULL, NULL, NULL, "--max",
     "turn off composition bias filter", 7},
    /* Other options */
    {"--nonull2", eslARG_NONE, NULL, NULL, NULL, NULL, NULL, NULL,
     "turn off biased composition score corrections", 12},
    {"-Z", eslARG_REAL, FALSE, NULL, "x>0", NULL, NULL, NULL,
     "set # of comparisons done, for E-value calculation", 12},
    {"--domZ", eslARG_REAL, FALSE, NULL, "x>0", NULL, NULL, NULL,
     "set # of significant seqs, for domain E-value calculation", 12},
    {"--seed", eslARG_INT, "42", NULL, "n>=0", NULL, NULL, NULL,
     "set RNG seed to <n> (if 0: one-time arbitrary seed)", 12},
    {"--qformat", eslARG_STRING, NULL, NULL, NULL, NULL, NULL, NULL,
     "assert input <seqfile> is in format <s>: no autodetection", 12},
    /* Not used, but retained because esl option-handling code errors if it isn't kept
       here.  Placed in group 99 so it doesn't print to help*/
    {"--notrans", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL,
     "don't show the translated DNA sequence in domain alignment",
     99}, /*for hmmscant */
    {"--vertcodon", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL,
     "show the DNA vertically in domain alignment", 99}, /*for hmmscant */
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
};

char usage[] = "[options] <hmmfile> <seqfile>";
char banner[] = "match DNA sequence against protein profile";

ESL_SQFILE *open_seqfile(const char *fpath, int fmt)
{
  ESL_SQFILE *seqfile = NULL;

  int status = esl_sqfile_Open(fpath, fmt, p7_SEQDBENV, &seqfile);
  if (status == eslENOTFOUND)
    p7_Fail("Failed to open sequence file %s for reading\n", fpath);
  else if (status == eslEFORMAT)
    p7_Fail("Sequence file %s is empty or misformatted\n", fpath);
  else if (status == eslEINVAL)
    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)
    p7_Fail("Unexpected error %d opening sequence file %s\n", status, fpath);

  return seqfile;
}

void workinfo_setup(struct WORKER_INFO *info, P7_HMMFILE *hfp, ESL_SQ *sq)
{
  info->th = p7_tophits_Create();
  info->pli = p7_pipeline_Create(info->cfg->go, 100, 100, FALSE, p7_SCAN_MODELS);
  info->pli->hfp = hfp;

  p7_pli_NewSeq(info->pli, sq);
  info->sq = sq;
}

void workinfo_reset(struct WORKER_INFO *info)
{
  p7_pipeline_Destroy(info->pli);
  p7_tophits_Destroy(info->th);
  info->pli = NULL;
  info->th = NULL;
}

void seq_fail(const char *fpath, int status, ESL_SQFILE *seqfile)
{
  if (status == eslEFORMAT)
    esl_fatal("Parse failed (sequence file %s):\n%s\n", fpath,
              esl_sqfile_GetErrorBuf(seqfile));
  else
    esl_fatal("Unexpected error %d reading sequence file %s", status, fpath);
}

P7_HMMFILE *open_hmmfile(const char *hfpath)
{
  P7_HMMFILE *hfp = NULL;
  char errbuf[eslERRBUFSIZE];
  int status = p7_hmmfile_OpenE(hfpath, p7_HMMDBENV, &hfp, errbuf);
  if (status == eslENOTFOUND)
    p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n",
            hfpath, errbuf);
  else if (status == eslEFORMAT)
    p7_Fail("File format problem, trying to open HMM file %s.\n%s\n", hfpath, errbuf);
  else if (status != eslOK)
    p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n", status, hfpath,
            errbuf);
  if (!hfp->is_pressed)
    p7_Fail("Failed to open binary auxfiles for %s: use hmmpress first\n", hfp->fname);
  return hfp;
}

void hmm_read_fail(const char *fpath, int status)
{
  switch (status) {
    case eslEOD:
      p7_Fail("read failed, HMM file %s may be truncated?", fpath);
      break;
    case eslEFORMAT:
      p7_Fail("bad file format in HMM file %s", fpath);
      break;
    case eslEINCOMPAT:
      p7_Fail("HMM file %s contains different alphabets", fpath);
      break;
    case eslEOF:
      p7_Fail("Unexpected end of HMM file %s", fpath);
      break;
    default:
      p7_Fail("Unexpected error (%d) in reading HMMs from %s", status, fpath);
  }
}

P7_HMM *read_hmm_model(const char *fpath)
{
  P7_HMM *hmm = NULL;
  ESL_ALPHABET *abc = NULL;
  P7_HMMFILE *hfp = open_hmmfile(fpath);
  int hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
  if (hstatus) hmm_read_fail(fpath, hstatus);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  return hmm;
}

void hmm_scan(struct WORKER_INFO *info, P7_HMMFILE *hfp)
{
  int status;

  P7_OPROFILE *om = NULL;
  ESL_ALPHABET *abc = NULL;
  P7_HMM *hmm = read_hmm_model(info->cfg->hfpath);

  while ((status = p7_oprofile_ReadMSV(hfp, &abc, &om)) == eslOK) {
    p7_pli_NewModel(info->pli, om, info->bg);
    p7_bg_SetLength(info->bg, info->sq->n);
    p7_oprofile_ReconfigLength(om, info->sq->n);

    P7_PROFILE *gm = p7_profile_Create(hmm->M, abc);
    p7_ProfileConfig(hmm, info->bg, gm, info->sq->n, p7_LOCAL);
    P7_GMX *gx = p7_gmx_Create(gm->M, info->sq->L);
    status = p7_PipelineHorta(info->pli, gm, om, gx, info->bg, info->sq, NULL, info->th,
                              NULL);
    p7_gmx_Destroy(gx);
    p7_profile_Destroy(gm);

    if (status == eslEINVAL) p7_Fail(info->pli->errbuf);

    p7_oprofile_Destroy(om);
    p7_pipeline_Reuse(info->pli);
  }

  esl_alphabet_Destroy(abc);
  p7_hmm_Destroy(hmm);

  if (status == eslEFORMAT)
    p7_Fail("bad file format in HMM file %s", info->cfg->hfpath);
  if (status == eslEINCOMPAT)
    p7_Fail("HMM file %s contains different alphabets", info->cfg->hfpath);
  if (status != eslOK && status != eslEOF)
    p7_Fail("Unexpected error in reading HMMs from %s", info->cfg->hfpath);
}

int show_seq_desc(FILE *ofp, const ESL_SQ *sq)
{
  int status = eslOK;
  FPRINTF(ofp, "Query:       %s  [L=%ld]\n", sq->name, (long)sq->n);
  if (sq->acc[0]) FPRINTF(ofp, "Accession:   %s\n", sq->acc);
  if (sq->desc[0]) FPRINTF(ofp, "Description: %s\n", sq->desc);
  return status;
}

int show_scan_results(FILE *ofp, struct WORKER_INFO *info, ESL_STOPWATCH *w)
{
  p7_tophits_SortBySortkey(info->th);
  p7_tophits_Threshold(info->th, info->pli);

  p7_tophits_Targets(ofp, info->th, info->pli, info->cfg->textw);
  FPUTS("\n\n", ofp);
  p7_tophits_Domains(ofp, info->th, info->pli, info->cfg->textw);
  FPUTS("\n\n", ofp);

  esl_stopwatch_Stop(w);
  p7_pli_Statistics(ofp, info->pli, w);
  FPUTS("//\n", ofp);
  fflush(ofp);
  return eslOK;
}

int foreach_seq(struct WORKER_INFO *info, ESL_SQFILE *sfp, ESL_SQ *sq)
{
  FILE *ofp = stdout;
  ESL_STOPWATCH *w = esl_stopwatch_Create();
  P7_HMMFILE *hfp = NULL;

  int sstatus = eslOK;
  int status = eslOK;
  int nquery = 0;

  while ((sstatus = esl_sqio_Read(sfp, sq)) == eslOK) {
    nquery++;
    esl_stopwatch_Start(w);

    if ((status = show_seq_desc(ofp, sq))) return status;

    hfp = open_hmmfile(info->cfg->hfpath);
    workinfo_setup(info, hfp, sq);

    hmm_scan(info, hfp);

    status = show_scan_results(ofp, info, w);
    if (status) return status;

    workinfo_reset(info);
    p7_hmmfile_Close(hfp);
    esl_sq_Reuse(sq);
  }
  if (sstatus != eslEOF) seq_fail(info->cfg->sfpath, sstatus, sfp);

  esl_stopwatch_Destroy(w);
  return status;
}

int get_seqfile_format(const ESL_GETOPTS *go)
{
  int seqfile_fmt = eslSQFILE_UNKNOWN;

  if (esl_opt_IsOn(go, "--qformat")) {
    seqfile_fmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--qformat"));
    if (seqfile_fmt == eslSQFILE_UNKNOWN)
      p7_Fail("%s is not a recognized input sequence file format\n",
              esl_opt_GetString(go, "--qformat"));
  }
  return seqfile_fmt;
}

ESL_ALPHABET *read_hmm_alphabet(const char *fpath)
{
  P7_HMM *hmm = NULL;
  ESL_ALPHABET *abc = NULL;
  P7_HMMFILE *hfp = open_hmmfile(fpath);
  int hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
  if (hstatus) hmm_read_fail(fpath, hstatus);
  p7_hmmfile_Close(hfp);
  p7_hmm_Destroy(hmm);
  return abc;
}

int show_hdr(FILE *ofp, const ESL_GETOPTS *go, const char *hfpath, const char *sfpath)
{
  p7_banner(ofp, go->argv[0], banner);
  const char *s = NULL;

  FPRINTF(ofp, "# Protein HMM file:         %s\n", hfpath);
  FPRINTF(ofp, "# DNA sequence file:        %s\n", sfpath);

  if (esl_opt_IsUsed(go, "--qformat"))
    FPRINTF(ofp, "# input seqfile format asserted:   %s\n",
            esl_opt_GetString(go, "--qformat"));

  s = "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n";
  FPUTS(s, ofp);

  return eslOK;
}

int serial_master(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int status = eslOK;
  show_hdr(stdout, go, cfg->hfpath, cfg->sfpath);

  ESL_ALPHABET *abc = read_hmm_alphabet(cfg->hfpath);
  ESL_SQFILE *sfp = open_seqfile(cfg->sfpath, get_seqfile_format(go));
  esl_sqfile_SetDigital(sfp, abc);
  ESL_SQ *sq = esl_sq_CreateDigital(abc);

  struct WORKER_INFO info = {
      .cfg = cfg,
      .bg = p7_bg_Create(abc),
      .pli = NULL,
      .th = NULL,
      .om = NULL,
      .sq = NULL,
  };

  foreach_seq(&info, sfp, sq);

  esl_sqfile_Close(sfp);
  esl_alphabet_Destroy(abc);
  return status;
}

void usage_err(const ESL_GETOPTS *go)
{
  int status = eslFAIL;

  esl_usage(stdout, go->argv[0], usage);
  XPUTS("\nwhere most common options are:");

  esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
  XPRINTF("\nTo see more help on available options, do %s -h.\n\n", go->argv[0]);

ERROR:
  exit(status);
}

void show_help(const ESL_GETOPTS *go)
{
  int status = eslOK;
  p7_banner(stdout, go->argv[0], banner);
  esl_usage(stdout, go->argv[0], usage);

  XPUTS("\nBasic options:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80);

  XPUTS("\nOptions directing output:");
  esl_opt_DisplayHelp(stdout, go, 2, 2, 80);

  XPUTS("\nOptions controlling reporting thresholds:");
  esl_opt_DisplayHelp(stdout, go, 4, 2, 80);

  XPUTS("\nOptions controlling inclusion (significance) thresholds:");
  esl_opt_DisplayHelp(stdout, go, 5, 2, 80);

  XPUTS("\nOptions controlling model-specific thresholding:");
  esl_opt_DisplayHelp(stdout, go, 6, 2, 80);

  XPUTS("\nOptions controlling acceleration heuristics:");
  esl_opt_DisplayHelp(stdout, go, 7, 2, 80);

  XPUTS("\nOther expert options:");
  esl_opt_DisplayHelp(stdout, go, 12, 2, 80);

  XPUTS("\nTranslation options:");
  esl_opt_DisplayHelp(stdout, go, 15, 2, 80);

ERROR:
  exit(status);
}

int cmdline(const ESL_GETOPTS *go, char **ret_hmm_fpath, char **ret_seq_fpath)
{
  if (esl_opt_GetBoolean(go, "-h") == TRUE) show_help(go);

  if (esl_opt_ArgNumber(go) != 2) {
    PUTS("Incorrect number of command line arguments.");
    usage_err(go);
  }

  if ((*ret_hmm_fpath = esl_opt_GetArg(go, 1)) == NULL) {
    PUTS("Failed to get <hmmfile> argument on command line.");
    usage_err(go);
  }

  if ((*ret_seq_fpath = esl_opt_GetArg(go, 2)) == NULL) {
    PUTS("Failed to get <seqfile> argument on command line.");
    usage_err(go);
  }

  return eslOK;
}

ESL_GETOPTS *create_getopts(int argc, char **argv)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  int status = eslOK;

  if (esl_opt_ProcessEnvironment(go)) {
    XPRINTF("Failed to process environment: %s\n", go->errbuf);
    usage_err(go);
  }

  if (esl_opt_ProcessCmdline(go, argc, argv)) {
    XPRINTF("Failed to parse command line: %s\n", go->errbuf);
    usage_err(go);
  }

  if (esl_opt_VerifyConfig(go)) {
    XPRINTF("Failed to parse command line: %s\n", go->errbuf);
    usage_err(go);
  }

  return go;

ERROR:
  exit(status);
}

int main(int argc, char **argv)
{
  impl_Init();
  p7_FLogsumInit();

  ESL_GETOPTS *go = create_getopts(argc, argv);

  struct cfg_s cfg = {.go = go, .sfpath = NULL, .hfpath = NULL, .textw = 0};

  int status = cmdline(go, &cfg.hfpath, &cfg.sfpath);
  if (status) return status;

  status = serial_master(go, &cfg);
  if (status) return status;

  esl_getopts_Destroy(go);
  return status;
}
