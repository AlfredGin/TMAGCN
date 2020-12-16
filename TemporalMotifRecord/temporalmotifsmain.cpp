// temporalmotifsmain.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include "temporalmotifs.h"
#include <stdio.h>
#include<string>
#include<sstream>
#include <vector>
#ifdef USE_OPENMP
#include <omp.h>
#endif
using namespace std;

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("Temporalmotifs. build: %s, %s. Time: %s",
			 __TIME__, __DATE__, TExeTm::GetCurTm()));  
  TExeTm ExeTm;  
  Try

  const TStr temporal_graph_filename =
    /*Env.GetIfArgPrefixStr("-i:", "example-temporal-graph.txt",
			  "input directed temporal graph file");*/
    /*Env.GetIfArgPrefixStr("-i:", "input-graph.txt",
              "Input directed temporal graph file");*/
      Env.GetIfArgPrefixStr("-i:", "socialevol-temporal-graph20081001.txt",
          "input directed temporal graph file");
      /*Env.GetIfArgPrefixStr("-i:", "Github-temporal-graph-sub.txt",
          "input directed temporal graph file");*/
      /*Env.GetIfArgPrefixStr("-i:", "CollegeMsg.txt",
              "input directed temporal graph file");*/
      /*Env.GetIfArgPrefixStr("-i:", "email-Eu-core-temporal-Dept3_sub.txt",
      "input directed temporal graph file");*/
  const TStr output = 
    /*Env.GetIfArgPrefixStr("-o:", "input-graph-motif-counts.txt",
			  "Output file in which to write counts");*/
      Env.GetIfArgPrefixStr("-o:", "socialevol-timewindow-1h-motif-counts.txt",
          "Output file in which to write counts");
      /*Env.GetIfArgPrefixStr("-o:", "Github-timewindow-1m-motif-counts.txt",
          "Output file in which to write counts");*/
    /*Env.GetIfArgPrefixStr("-o:", "email-Eu-core-temporal-Dept3_sub-timewindow-1d-motif-counts.txt",
              "Output file in which to write counts");*/
  /*Env.GetIfArgPrefixStr("-o:", "CollegeMsg-timewindow-10min-motif-counts.txt",
      "Output file in which to write counts");*/
  const TStr outputRCount =
      /*Env.GetIfArgPrefixStr("-o:", "input-graph-motif-counts.txt",
                "Output file in which to write counts");*/
      Env.GetIfArgPrefixStr("-o:", "socialevol-timewindow-1h-motif-Rcounts.txt",
          "Output file in which to write counts");
      /*Env.GetIfArgPrefixStr("-o:", "Github-timewindow-1m-motif-Rcounts.txt",
          "Output file in which to write counts");*/
      /*Env.GetIfArgPrefixStr("-o:", "email-Eu-core-temporal-Dept3_sub-timewindow-1d-motif-Rcounts.txt",
          "Output file in which to write counts");*/
      /*Env.GetIfArgPrefixStr("-o:", "CollegeMsg-timewindow-10min-motif-Rcounts.txt",
          "Output file in which to write counts");*/
  const TStr outputR =
      /*Env.GetIfArgPrefixStr("-or:", "input-graph-motif-records.txt",
          "Output file in which to write motif records");*/
      Env.GetIfArgPrefixStr("-or:", "socialevol-timewindow-1h-motif-records.txt",
          "Output file in which to write motif records");
      /*Env.GetIfArgPrefixStr("-or:", "Github-timewindow-1m-motif-records.txt",
          "Output file in which to write motif records");*/
      /*Env.GetIfArgPrefixStr("-or:", "email-Eu-core-temporal-Dept3_sub-timewindow-1d-motif-records.txt",
          "Output file in which to write motif records");*/
  /*Env.GetIfArgPrefixStr("-or:", "CollegeMsg-timewindow-10min-motif-records.txt",
      "Output file in which to write motif records");*/


  const TFlt delta =
    Env.GetIfArgPrefixFlt("-delta:", 60*60, "Time window delta");
  const int num_threads =
    Env.GetIfArgPrefixInt("-nt:", 8, "Number of threads for parallelization");

  /*TVec<TStr> motifRecordFile;
  TStr motifRecordFileName;
  for (TInt i = 0; i < 5; i++){
      for (TInt j = 0; j < 5; j++) {
          motifRecordFileName = "MotifRecord" + i.GetStr()+ j.GetStr();
          const TStr motifRecordFileName =
              Env.GetIfArgPrefixStr("-or:", motifRecordFileName,
                  "Output file in which to write motif records");
          motifRecordFile.Add(motifRecordFileName);
      }
  }*/

#ifdef USE_OPENMP
  omp_set_num_threads(num_threads);
#endif

  // Count all 2-node and 3-node temporal motifs with 3 temporal edges
  TempMotifCounter tmc(temporal_graph_filename);
  Counter2D counts;
  //TVec<TStr> motifRecordFile;
  TVec<TVec<TStr>> motifRecordFile(44);
  TStr motifRecordIndex;
  /*for (TInt i = 0; i < 6; i++) {
      for (TInt j = 0; j < 6; j++) {
          motifRecordIndex = "MotifRecord" + i.GetStr() + j.GetStr();
          TVec<TStr> motifRecordIndex;
          motifRecordFile.Add(motifRecordIndex);
      }
  }*/

  tmc.Count3TEdge23Node(delta, counts, &motifRecordFile);
  FILE* output_file = fopen(output.CStr(), "wt");
  for (int i = 0; i < counts.m(); i++) {
    for (int j = 0; j < counts.n(); j++) {
      int count = counts(i, j);
      fprintf(output_file, "%d", count);
      if (j < counts.n() - 1) { fprintf(output_file, " "); }
    }
    fprintf(output_file, "\n");
  }

  FILE* output_fileRC = fopen(outputRCount.CStr(), "wt");
  for (int i = 0; i < counts.m(); i++) {
      for (int j = 0; j < counts.n(); j++) {
          fprintf(output_fileRC, "%d", motifRecordFile[i*6+j].Len());
          if (j < counts.n() - 1) { fprintf(output_fileRC, " "); }
      }
      fprintf(output_fileRC, "\n");
  }

  FILE* output_fileR = fopen(outputR.CStr(), "wt");
  //printf("\nmotifRecordFile.Len():%d", motifRecordFile.Len());
  for (int i = 0; i < motifRecordFile.Len(); i++) {
      printf("\nmotifRecordFile[%d].Len():%d", i,motifRecordFile[i].Len());
      for (int j = 0; j < motifRecordFile[i].Len(); j++) {
          fprintf(output_fileR, "%s", motifRecordFile[i][j].CStr());
      }
  }

  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(),
	 TSecTm::GetCurTm().GetTmStr().CStr());
  return 0;
}
