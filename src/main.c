/*
OpenPEEC Version 1.0.0

準静的 PEEC (部分要素等価回路) 回路ソルバー : solver main
*/

#include "peec.h"

#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static void args(int, char *[], int *, char []);
static void monitor(FILE *fp, const char *str);
static double walltime(void);

int main(int argc, char *argv[])
{
	const char errfmt[] = "*** file %s open error.\n";
	char str[BUFSIZ];
	FILE *fp_in = NULL, *fp_log = NULL;
	peec_t peec;

	// arguments
	int nthread = 1;
	char fn_in[BUFSIZ] = "";
	args(argc, argv, &nthread, fn_in);

	// set number of threads
#ifdef _OPENMP
	omp_set_num_threads(nthread);
#endif

	const double cpu0 = walltime();

	// input data
	if ((fp_in = fopen(fn_in, "r")) == NULL) {
		printf(errfmt, fn_in);
		exit(1);
	}

	// open log file
	if ((fp_log = fopen(FN_LOG, "w")) == NULL) {
		printf(errfmt, FN_LOG);
		exit(1);
	}

	// logo
	sprintf(str, "<<< %s Ver.%d.%d.%d >>>", PROGRAM, VERSION_MAJOR, VERSION_MINOR, VERSION_BUILD);
	monitor(fp_log, str);
#ifdef _OPENMP
	sprintf(str, "CPU, thread=%d", nthread);
#else
	sprintf(str, "CPU, no OpenMP");
#endif
	monitor(fp_log, str);

	// input
	if (input_data(fp_in, &peec)) {
		fprintf(fp_log, "%s\n", "*** input data error");
		fclose(fp_log);
		fclose(fp_in);
		exit(1);
	}
	fclose(fp_in);
	printf("title = %s\n", peec.title);
	fprintf(fp_log, "title = %s\n", peec.title);
	snprintf(str, sizeof(str), "R=%d C=%d L=%d K=%d source=%d wire=%d port=%d frequency=%d",
		peec.nres, peec.ncap, peec.nind, peec.nmut, peec.nsrc, peec.nwire, peec.nport, peec.nfreq);
	monitor(fp_log, str);

	// setup : ワイヤ分割 -> 部分インダクタンス -> 電位係数 -> MNA 番号付け
	if (wire_build(&peec, fp_log)) {
		fclose(fp_log);
		exit(1);
	}
	lp_fill(&peec, fp_log);
	if (pot_fill(&peec, fp_log)) {
		fclose(fp_log);
		exit(1);
	}
	if (mna_numbering(&peec, fp_log)) {
		fclose(fp_log);
		exit(1);
	}

	const double cpu1 = walltime();

	// solve
	if (solve(&peec, fp_log)) {
		fclose(fp_log);
		exit(1);
	}

	const double cpu2 = walltime();

	// output
	output_zin(&peec, fp_log);
	if (output_csv(&peec, FN_CSV)) {
		fclose(fp_log);
		exit(1);
	}
	sprintf(str, "output filename : %s, %s", FN_LOG, FN_CSV);
	monitor(fp_log, str);

	// cpu time
	sprintf(str, "cpu time : setup %.3f s, solve %.3f s", cpu1 - cpu0, cpu2 - cpu1);
	monitor(fp_log, str);

	monitor(fp_log, "=== normal end ===");

	fclose(fp_log);

	return 0;
}


static void args(int argc, char *argv[], int *nthread, char fn_in[])
{
	const char usage[] = "Usage : peec [-n <thread>] <datafile>";

	if (argc < 2) {
		printf("%s\n", usage);
		exit(0);
	}

	while (--argc) {
		++argv;
		if (!strcmp(*argv, "-n")) {
			if (--argc) {
				*nthread = atoi(*++argv);
				if (*nthread < 1) *nthread = 1;
			}
			else {
				break;
			}
		}
		else if (!strcmp(*argv, "--help")) {
			printf("%s\n", usage);
			exit(0);
		}
		else if (!strcmp(*argv, "--version")) {
			printf("%s Ver.%d.%d.%d\n", PROGRAM, VERSION_MAJOR, VERSION_MINOR, VERSION_BUILD);
			exit(0);
		}
		else {
			strcpy(fn_in, *argv);
		}
	}

	if (fn_in[0] == '\0') {
		printf("%s\n", usage);
		exit(0);
	}
}


// stdout とログの両方に出力
static void monitor(FILE *fp, const char *str)
{
	printf("%s\n", str);
	fflush(stdout);
	fprintf(fp, "%s\n", str);
	fflush(fp);
}


static double walltime(void)
{
#ifdef _OPENMP
	return omp_get_wtime();
#else
	return (double)clock() / CLOCKS_PER_SEC;
#endif
}
