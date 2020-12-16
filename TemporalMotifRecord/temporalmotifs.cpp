#include <iostream>
#include <algorithm>
#include "Snap.h"
#include "temporalmotifs.h"
using namespace std;
///////////////////////////////////////////////////////////////////////////////
// Initialization and helper methods for TempMotifCounter
TempMotifCounter::TempMotifCounter(const TStr& filename) {
  // First load the static graph
  static_graph_ = TSnap::LoadEdgeList<PNGraph>(filename, 0, 1);
  int max_nodes = static_graph_->GetMxNId();
  cout << "maxnodes:"<<max_nodes<< endl;
  temporal_data_ = TVec< THash<TInt, TIntV> >(max_nodes);

  // Formulate input File Format:
  //   source_node destination_node timestamp
  TTableContext context;
  Schema temp_graph_schema;
  temp_graph_schema.Add(TPair<TStr,TAttrType>("source", atInt));
  temp_graph_schema.Add(TPair<TStr,TAttrType>("destination", atInt));
  temp_graph_schema.Add(TPair<TStr,TAttrType>("time", atInt));

  // Load the temporal graph
  PTable data_ptr = TTable::LoadSS(temp_graph_schema, filename, &context, ' ');
  TInt src_idx = data_ptr->GetColIdx("source");
  TInt dst_idx = data_ptr->GetColIdx("destination");
  TInt tim_idx = data_ptr->GetColIdx("time");
  for (TRowIterator RI = data_ptr->BegRI(); RI < data_ptr->EndRI(); RI++) {
    TInt row_idx = RI.GetRowIdx();
    int src = data_ptr->GetIntValAtRowIdx(src_idx, row_idx).Val;
    int dst = data_ptr->GetIntValAtRowIdx(dst_idx, row_idx).Val;
    int tim = data_ptr->GetIntValAtRowIdx(tim_idx, row_idx).Val;
    // Do not include self loops as they do not appear in the definition of
    // temporal motifs.
    if (src != dst) { temporal_data_[src](dst).Add(tim); }
  }
}

void TempMotifCounter::GetAllNodes(TIntV& nodes) {
  nodes = TIntV();
  for (TNGraph::TNodeI it = static_graph_->BegNI();
       it < static_graph_->EndNI(); it++) {
    nodes.Add(it.GetId());
  }
}

bool TempMotifCounter::HasEdges(int u, int v) {
  return temporal_data_[u].IsKey(v);
}

void TempMotifCounter::GetAllNeighbors(int node, TIntV& nbrs) {
  nbrs = TIntV();
  TNGraph::TNodeI NI = static_graph_->GetNI(node);
  for (int i = 0; i < NI.GetOutDeg(); i++) { nbrs.Add(NI.GetOutNId(i)); }
  for (int i = 0; i < NI.GetInDeg(); i++) {
    int nbr = NI.GetInNId(i);
    if (!NI.IsOutNId(nbr)) { nbrs.Add(nbr); }
  }
}

void TempMotifCounter::GetAllStaticTriangles(TIntV& Us, TIntV& Vs, TIntV& Ws) {
  Us.Clr();
  Vs.Clr();
  Ws.Clr();
  // Get degree ordering of the graph
  int max_nodes = static_graph_->GetMxNId();
  TVec<TIntPair> degrees(max_nodes);
  degrees.PutAll(TIntPair(0, 0));
  // Set the degree of a node to be the number of nodes adjacent to the node in
  // the undirected graph.
  TIntV nodes;
  GetAllNodes(nodes);
  #pragma omp parallel for schedule(dynamic)  
  for (int node_id = 0; node_id < nodes.Len(); node_id++) {
    int src = nodes[node_id];
    TIntV nbrs;
    GetAllNeighbors(src, nbrs);
    degrees[src] = TIntPair(nbrs.Len(), src);
  }
  degrees.Sort();
  TIntV order = TIntV(max_nodes);
  #pragma omp parallel for schedule(dynamic)  
  for (int i = 0; i < order.Len(); i++) {
    order[degrees[i].Dat] = i;
  }

  // Get triangles centered at a given node where that node is the smallest in
  // the degree ordering.
  #pragma omp parallel for schedule(dynamic)  
  for (int node_id = 0; node_id < nodes.Len(); node_id++) {
    int src = nodes[node_id];
    int src_pos = order[src];
    
    // Get all neighbors who come later in the ordering
    TIntV nbrs;
    GetAllNeighbors(src, nbrs);    
    TIntV neighbors_higher;
    for (int i = 0; i < nbrs.Len(); i++) {
      int nbr = nbrs[i];
      if (order[nbr] > src_pos) { neighbors_higher.Add(nbr); }
    }

    for (int ind1 = 0; ind1 < neighbors_higher.Len(); ind1++) {
      for (int ind2 = ind1 + 1; ind2 < neighbors_higher.Len(); ind2++) {
        int dst1 = neighbors_higher[ind1];
        int dst2 = neighbors_higher[ind2];
        // Check for triangle formation
        if (static_graph_->IsEdge(dst1, dst2) || static_graph_->IsEdge(dst2, dst1)) {
          #pragma omp critical
          {
            Us.Add(src);
            Vs.Add(dst1);
            Ws.Add(dst2);
          }
        }
      }
    }
  }
}

void TempMotifCounter::Count3TEdge23Node(double delta, Counter2D& counts, TVec<TVec<TStr>>* motifRecordFile) {
  // This is imply a wrapper function around the counting methods to produce
  // counts in the same way that they were represented in the paper.  This makes
  // it easy to reproduce results and allow SNAP users to make the same
  // measurements on their temporal network data.
  counts = Counter2D(6, 6);
  
  Counter2D edge_counts;
  /*TIntV a;
  for (int i = 0; i < 5; i++) {
      a.Add(TInt(i));
      printf("\na: %d", i);
  }
  TIntV b;
  for (int i = 3; i < 7; i++) {
      b.Add(TInt(i));
      printf("\nb: %d", i);
  }
  a.Diff(b);
  for (int i = 0; i < a.Len(); i++) {
      printf("\naDiff: %d", a[i]);
  }*/


  Count3TEdge2Node(delta, edge_counts, motifRecordFile);
  counts(4, 0) = edge_counts(0, 0);
  counts(4, 1) = edge_counts(0, 1);
  counts(5, 0) = edge_counts(1, 0);
  counts(5, 1) = edge_counts(1, 1);

  Counter3D pre_counts, pos_counts, mid_counts;
  //Count3TEdge3NodeStars(delta, pre_counts, pos_counts, mid_counts, motifRecordFile);
  Count3TEdge3NodeStarsNaive(delta, pre_counts, pos_counts, mid_counts, motifRecordFile);
  counts(0, 0) = mid_counts(1, 1, 1);
  counts(0, 1) = mid_counts(1, 1, 0);
  counts(0, 4) = pos_counts(1, 1, 0);
  counts(0, 5) = pos_counts(1, 1, 1);
  counts(1, 0) = mid_counts(1, 0, 1);
  counts(1, 1) = mid_counts(1, 0, 0);
  counts(1, 4) = pos_counts(1, 0, 0);
  counts(1, 5) = pos_counts(1, 0, 1);
  counts(2, 0) = mid_counts(0, 1, 0);
  counts(2, 1) = mid_counts(0, 1, 1);
  counts(2, 2) = pos_counts(0, 1, 0);
  counts(2, 3) = pos_counts(0, 1, 1);
  counts(3, 0) = mid_counts(0, 0, 0);
  counts(3, 1) = mid_counts(0, 0, 1);
  counts(3, 2) = pos_counts(0, 0, 0);
  counts(3, 3) = pos_counts(0, 0, 1);
  counts(4, 2) = pre_counts(0, 1, 0);
  counts(4, 3) = pre_counts(0, 1, 1);
  counts(4, 4) = pre_counts(1, 0, 0);
  counts(4, 5) = pre_counts(1, 0, 1);
  counts(5, 2) = pre_counts(0, 0, 0);
  counts(5, 3) = pre_counts(0, 0, 1);
  counts(5, 4) = pre_counts(1, 1, 0);
  counts(5, 5) = pre_counts(1, 1, 1);  

  Counter3D triad_counts;
  //Count3TEdgeTriads(delta, triad_counts);
  Count3TEdgeTriadsNaive(delta, triad_counts, motifRecordFile);
  counts(0, 2) = triad_counts(0, 0, 0);
  counts(0, 3) = triad_counts(0, 0, 1);
  counts(1, 2) = triad_counts(0, 1, 0);
  counts(1, 3) = triad_counts(0, 1, 1);
  counts(2, 4) = triad_counts(1, 0, 0);
  counts(2, 5) = triad_counts(1, 0, 1);
  counts(3, 4) = triad_counts(1, 1, 0);
  counts(3, 5) = triad_counts(1, 1, 1);
}

///////////////////////////////////////////////////////////////////////////////
// Two-node (static edge) counting methods
void TempMotifCounter::Count3TEdge2Node(double delta, Counter2D& counts, TVec<TVec<TStr>>* motifRecordFile) {
  // Get a vector of undirected edges (so we can use openmp parallel for over it)
  TVec<TIntPair> undir_edges;
  
  for (TNGraph::TEdgeI it = static_graph_->BegEI(); it < static_graph_->EndEI(); it++) {
      int src = it.GetSrcNId();
      int dst = it.GetDstNId();
      // Only consider undirected edges
      if (src < dst || (dst < src && !static_graph_->IsEdge(dst, src))) {
          undir_edges.Add(TIntPair(src, dst));
      }
  }
  printf("\nundir_edges: %d \n", undir_edges.Len());
  counts = Counter2D(2, 2);
  #pragma omp parallel for schedule(dynamic)

  for (int i = 0; i < undir_edges.Len(); i++) {
    TIntPair edge = undir_edges[i];
    //printf("\nundir_edges: %d+%d", edge.Key, edge.Dat);
    Counter3D local;
    Count3TEdge2Node(edge.Key, edge.Dat, delta, local, motifRecordFile);//s:1+d:2  s:1+d:3  s:1+d:4  s:3+d:4  
    #pragma omp critical
    {
      counts(0, 0) += local(0, 1, 0) + local(1, 0, 1);  // M_{5,1}
      counts(0, 1) += local(1, 0, 0) + local(0, 1, 1);  // M_{5,2}
      counts(1, 0) += local(0, 0, 0) + local(1, 1, 1);  // M_{6,1}
      counts(1, 1) += local(0, 0, 1) + local(1, 1, 0);  // M_{6,2}
    }
  }
}

void TempMotifCounter::Count3TEdge2Node(int u, int v, double delta,
                                        Counter3D& counts, TVec<TVec<TStr>>* motifRecordFile) {
  // Sort event list by time
  //TVec<TIntPair> combined;
  //TThree<TIntV:in_out,TIntV:timestamps,TIntPair:uv>
  TVec<TThree> combined;
  AddStarEdges(combined, u, v, 0);
  AddStarEdges(combined, v, u, 1);
  combined.Sort();

  // Get the counts
  ThreeTEdgeMotifCounter counter(2);
  TIntV in_out(combined.Len());
  TIntV timestamps(combined.Len());
  TVec<TIntPair> uv(combined.Len());
  //MotifR=TQuad<TInt, TInt, TIntPair, TInt> TRecord=<lable, event, TIntPair uv, TInt timestamps>
  //TVec<TRecord> MotifR;
  for (int k = 0; k < combined.Len(); k++) {
    timestamps[k] = combined[k].Val1;
    in_out[k] = combined[k].Val2;
    uv[k] = combined[k].Val3;
  }
  counter.Count(in_out, timestamps, uv, delta, counts, motifRecordFile);
  //printf("\ncombined: %d \n", combined.Len());
}

///////////////////////////////////////////////////////////////////////////////
// Star counting methods
//void TempMotifCounter::AddStarEdges(TVec<TIntPair>& combined, int u, int v,
//                                    int key) {
//  if (HasEdges(u, v)) {
//    const TIntV& timestamps = temporal_data_[u].GetDat(v);
//    for (int i = 0; i < timestamps.Len(); i++) {
//      combined.Add(TIntPair(timestamps[i], key));
//    }
//    printf("\norigin_timestamps: %d \n", timestamps.Len());  //s:1+d:2 len=1 s:1+d:3 len=4 s:1+d:4 len=1 s:3+d:4 len=1
//  }
//}

void TempMotifCounter::AddStarEdges(TVec<TThree>& combined, int u, int v,
    int key) {
    if (HasEdges(u, v)) {
        const TIntV& timestamps = temporal_data_[u].GetDat(v);
        for (int i = 0; i < timestamps.Len(); i++) {
            combined.Add(TThree(timestamps[i], key, TIntPair(u, v)));
        }
        //printf("\norigin_timestamps: %d \n", timestamps.Len());  //s:1+d:2 len=1 s:1+d:3 len=4 s:1+d:4 len=1 s:3+d:4 len=1
    }
}

void TempMotifCounter::Count3TEdge3NodeStarsNaive(
        double delta, Counter3D& pre_counts, Counter3D& pos_counts,
        Counter3D& mid_counts, TVec<TVec<TStr>>* motifRecordFile) {
  TIntV centers;
  GetAllNodes(centers);
  pre_counts = Counter3D(2, 2, 2);
  pos_counts = Counter3D(2, 2, 2);
  mid_counts = Counter3D(2, 2, 2);
  // Get counts for each node as the center
  #pragma omp parallel for schedule(dynamic)
  for (int c = 0; c < centers.Len(); c++) {
    // Gather all adjacent events
    int center = centers[c];
    TIntV nbrs;
    //MotifR=TQuad<TInt, TInt, TIntPair, TInt> TRecord=<lable, event, TIntPair uv, TInt timestamps>
    TVec<TRecord> MotifR;
    GetAllNeighbors(center, nbrs);
    for (int i = 0; i < nbrs.Len(); i++) {
      for (int j = i + 1; j < nbrs.Len(); j++) {
        int nbr1 = nbrs[i];
        int nbr2 = nbrs[j];
        //TVec<TIntPair> combined;
        //TThree<TIntV:in_out,TIntV:timestamps,TIntPair:uv>
        TVec<TThree> combined;
        AddStarEdges(combined, center, nbr1, 0);
        AddStarEdges(combined, nbr1, center, 1);
        AddStarEdges(combined, center, nbr2, 2);
        AddStarEdges(combined, nbr2, center, 3);
        combined.Sort();
        ThreeTEdgeMotifCounter counter(4);
        TIntV edge_id(combined.Len());
        TIntV timestamps(combined.Len());
        TVec<TIntPair> uv(combined.Len());
        for (int k = 0; k < combined.Len(); k++) {
          timestamps[k] = combined[k].Val1;
          edge_id[k] = combined[k].Val2;
          uv[k] = combined[k].Val3;
        }
        Counter3D local;
        //counter.Count(edge_id, timestamps, delta, local);
        counter.Count(edge_id, timestamps, uv, delta, local, motifRecordFile);

        #pragma omp critical
        {  // Update with local counts
          for (int dir1 = 0; dir1 < 2; ++dir1) {
            for (int dir2 = 0; dir2 < 2; ++dir2) {
              for (int dir3 = 0; dir3 < 2; ++dir3) {
                pre_counts(dir1, dir2, dir3) +=
                  local(dir1, dir2, dir3 + 2) + local(dir1 + 2, dir2 + 2, dir3);
                pos_counts(dir1, dir2, dir3) +=
                  local(dir1, dir2 + 2, dir3 + 2) + local(dir1 + 2, dir2, dir3);
                mid_counts(dir1, dir2, dir3) +=
                  local(dir1, dir2 + 2, dir3) + local(dir1 + 2, dir2, dir3 + 2);
              }
            }
          }
        }
      }
    }
  }
}

void TempMotifCounter::AddStarEdgeData(TVec<TIntPair>& ts_indices,
                                       TVec<TThree>& events,
                                       int& index, int u, int v, int nbr, int key) {
  if (HasEdges(u, v)) {
    const TIntV& ts_vec = temporal_data_[u].GetDat(v);
    for (int j = 0; j < ts_vec.Len(); ++j) {
      ts_indices.Add(TIntPair(ts_vec[j], index));
      events.Add(TThree(nbr, key, TIntPair(u, v)));
      index++;
    }
  }
}

void TempMotifCounter::Count3TEdge3NodeStars(double delta, Counter3D& pre_counts,
                                             Counter3D& pos_counts,
                                             Counter3D& mid_counts, TVec<TVec<TStr>>* motifRecordFile) {
  TIntV centers;
  GetAllNodes(centers);
  pre_counts = Counter3D(2, 2, 2);
  pos_counts = Counter3D(2, 2, 2);
  mid_counts = Counter3D(2, 2, 2);
  // Get counts for each node as the center
  #pragma omp parallel for schedule(dynamic)  
  for (int c = 0; c < centers.Len(); c++) {
    // Gather all adjacent events
    int center = centers[c];
    TVec<TIntPair> ts_indices;
    TVec<TThree> events;
    TNGraph::TNodeI NI = static_graph_->GetNI(center);
    int index = 0;
    TIntV nbrs;
    GetAllNeighbors(center, nbrs);
    int nbr_index = 0;
    for (int i = 0; i < nbrs.Len(); i++) {
      int nbr = nbrs[i];
      AddStarEdgeData(ts_indices, events, index, center, nbr, nbr_index, 0);
      AddStarEdgeData(ts_indices, events, index, nbr, center, nbr_index, 1);
      nbr_index++;
    }
    ts_indices.Sort();
    // ts_indices:(Key:timestamps,Dat:index)  events:(Key:nbr_index,Dat:key(0:out,1:in))
    TIntV timestamps;
    TVec<TThree> ordered_events;
    for (int j = 0; j < ts_indices.Len(); j++) {
      timestamps.Add(ts_indices[j].Key);
      ordered_events.Add(events[ts_indices[j].Dat]);
    }
    
    ThreeTEdgeStarCounter tesc(nbr_index);
    // dirs: outgoing --> 0, incoming --> 1
    tesc.Count(ordered_events, timestamps, delta, motifRecordFile);
    #pragma omp critical
    { // Update counts
      for (int dir1 = 0; dir1 < 2; ++dir1) {
        for (int dir2 = 0; dir2 < 2; ++dir2) {
          for (int dir3 = 0; dir3 < 2; ++dir3) {
            pre_counts(dir1, dir2, dir3) += tesc.PreCount(dir1, dir2, dir3);
            pos_counts(dir1, dir2, dir3) += tesc.PosCount(dir1, dir2, dir3);
            mid_counts(dir1, dir2, dir3) += tesc.MidCount(dir1, dir2, dir3);
          }
        }
      }
    }

    // Subtract off edge-wise counts
    for (int nbr_id = 0; nbr_id < nbrs.Len(); nbr_id++) {
      int nbr = nbrs[nbr_id];
      Counter3D edge_counts;
      Count3TEdge2Node(center, nbr, delta, edge_counts, motifRecordFile);
      #pragma omp critical
      {
        for (int dir1 = 0; dir1 < 2; ++dir1) {
          for (int dir2 = 0; dir2 < 2; ++dir2) {
            for (int dir3 = 0; dir3 < 2; ++dir3) {
              pre_counts(dir1, dir2, dir3) -= edge_counts(dir1, dir2, dir3);
              pos_counts(dir1, dir2, dir3) -= edge_counts(dir1, dir2, dir3);
              mid_counts(dir1, dir2, dir3) -= edge_counts(dir1, dir2, dir3);
            }
          }
        }
      }
    }
  }
  (*motifRecordFile)[3 * 6 + 0].Diff((*motifRecordFile)[36 + 2 * 2 * 0 + 2 * 0 + 1 * 0]);
  (*motifRecordFile)[5 * 6 + 2].Diff((*motifRecordFile)[36 + 2 * 2 * 0 + 2 * 0 + 1 * 0]);
  (*motifRecordFile)[3 * 6 + 2].Diff((*motifRecordFile)[36 + 2 * 2 * 0 + 2 * 0 + 1 * 0]);

  (*motifRecordFile)[3 * 6 + 1].Diff((*motifRecordFile)[36 + 2 * 2 * 0 + 2 * 0 + 1 * 1]);
  (*motifRecordFile)[5 * 6 + 3].Diff((*motifRecordFile)[36 + 2 * 2 * 0 + 2 * 0 + 1 * 1]);
  (*motifRecordFile)[3 * 6 + 3].Diff((*motifRecordFile)[36 + 2 * 2 * 0 + 2 * 0 + 1 * 1]);

  (*motifRecordFile)[2 * 6 + 0].Diff((*motifRecordFile)[36 + 2 * 2 * 0 + 2 * 1 + 1 * 0]);
  (*motifRecordFile)[4 * 6 + 2].Diff((*motifRecordFile)[36 + 2 * 2 * 0 + 2 * 1 + 1 * 0]);
  (*motifRecordFile)[2 * 6 + 2].Diff((*motifRecordFile)[36 + 2 * 2 * 0 + 2 * 1 + 1 * 0]);

  (*motifRecordFile)[2 * 6 + 1].Diff((*motifRecordFile)[36 + 2 * 2 * 0 + 2 * 1 + 1 * 1]);
  (*motifRecordFile)[4 * 6 + 3].Diff((*motifRecordFile)[36 + 2 * 2 * 0 + 2 * 1 + 1 * 1]);
  (*motifRecordFile)[2 * 6 + 3].Diff((*motifRecordFile)[36 + 2 * 2 * 0 + 2 * 1 + 1 * 1]);

  (*motifRecordFile)[1 * 6 + 1].Diff((*motifRecordFile)[36 + 2 * 2 * 1 + 2 * 0 + 1 * 0]);
  (*motifRecordFile)[4 * 6 + 4].Diff((*motifRecordFile)[36 + 2 * 2 * 1 + 2 * 0 + 1 * 0]);
  (*motifRecordFile)[1 * 6 + 4].Diff((*motifRecordFile)[36 + 2 * 2 * 1 + 2 * 0 + 1 * 0]);

  (*motifRecordFile)[1 * 6 + 0].Diff((*motifRecordFile)[36 + 2 * 2 * 1 + 2 * 0 + 1 * 1]);
  (*motifRecordFile)[4 * 6 + 5].Diff((*motifRecordFile)[36 + 2 * 2 * 1 + 2 * 0 + 1 * 1]);
  (*motifRecordFile)[1 * 6 + 5].Diff((*motifRecordFile)[36 + 2 * 2 * 1 + 2 * 0 + 1 * 1]);

  (*motifRecordFile)[0 * 6 + 1].Diff((*motifRecordFile)[36 + 2 * 2 * 1 + 2 * 1 + 1 * 0]);
  (*motifRecordFile)[5 * 6 + 4].Diff((*motifRecordFile)[36 + 2 * 2 * 1 + 2 * 1 + 1 * 0]);
  (*motifRecordFile)[0 * 6 + 4].Diff((*motifRecordFile)[36 + 2 * 2 * 1 + 2 * 1 + 1 * 0]);

  (*motifRecordFile)[0 * 6 + 0].Diff((*motifRecordFile)[36 + 2 * 2 * 1 + 2 * 1 + 1 * 1]);
  (*motifRecordFile)[5 * 6 + 5].Diff((*motifRecordFile)[36 + 2 * 2 * 1 + 2 * 1 + 1 * 1]);
  (*motifRecordFile)[0 * 6 + 5].Diff((*motifRecordFile)[36 + 2 * 2 * 1 + 2 * 1 + 1 * 1]);
}

/////////////////////////////////////////////////////////////////////////////
 //Triad counting methods
void TempMotifCounter::Count3TEdgeTriadsNaive(double delta, Counter3D& counts, TVec<TVec<TStr>>* motifRecordFile) {
  TIntV Us, Vs, Ws;
  GetAllStaticTriangles(Us, Vs, Ws);
  counts = Counter3D(2, 2, 2);
  #pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < Us.Len(); i++) {
    int u = Us[i];
    int v = Vs[i];
    int w = Ws[i];
    // Gather all edges in triangle (u, v, w)
    int uv = 0, vu = 1, uw = 2, wu = 3, vw = 4, wv = 5;
    //TVec<TIntPair> combined;
    TVec<TThree> combined;
    AddStarEdges(combined, u, v, uv);
    AddStarEdges(combined, v, u, vu);
    AddStarEdges(combined, u, w, uw);
    AddStarEdges(combined, w, u, wu);
    AddStarEdges(combined, v, w, vw);
    AddStarEdges(combined, w, v, wv);        
    // Get the counts for this triangle
    combined.Sort();
    ThreeTEdgeMotifCounter counter(6);
    TIntV edge_id(combined.Len());
    TIntV timestamps(combined.Len());
    TVec<TIntPair> uvpair(combined.Len());
    for (int k = 0; k < combined.Len(); k++) {
      /*edge_id[k] = combined[k].Dat;
      timestamps[k] = combined[k].Key;*/
      edge_id[k] = combined[k].Val2;
      timestamps[k] = combined[k].Val1;
      uvpair[k] = combined[k].Val3;
    }
    Counter3D local;
    //counter.Count(edge_id, timestamps, delta, local);
    counter.Count(edge_id, timestamps, uvpair, delta, local, motifRecordFile);

    // Update the global counter with the various symmetries
    #pragma omp critical
    {
      // i --> j, k --> j, i --> k
      counts(0, 0, 0) += local(uv, wv, uw) + local(vu, wu, vw) + local(uw, vw, uv)
        + local(wu, vu, wv) + local(vw, uw, vu) + local(wv, uv, wu);
      // i --> j, k --> j, k --> i
      counts(0, 0, 1) += local(uv, wv, wu) + local(vu, wu, wv) + local(uw, vw, vu)
        + local(wu, vu, vw) + local(vw, uw, uv) + local(wv, uv, uw);
      // i --> j, j --> k, i --> k
      counts(0, 1, 0) += local(uv, vw, uw) + local(vu, uw, vw) + local(uw, wv, uv)
        + local(wu, uv, wv) + local(vw, wu, vu) + local(wv, vu, wu);
      // i --> j, j --> k, k --> i
      counts(0, 1, 1) += local(uv, vw, wu) + local(vu, uw, wv) + local(uw, wv, vu)
        + local(wu, uv, vw) + local(vw, wu, uv) + local(wv, vu, uw);
      // i --> j, k --> i, j --> k
      counts(1, 0, 0) += local(uv, wu, vw) + local(vu, wv, uw) + local(uw, vu, wv)
        + local(wu, vw, uv) + local(vw, uv, wu) + local(wv, uw, vu);
      // i --> j, k --> i, k --> j
      counts(1, 0, 1) += local(uv, wu, wv) + local(vu, wv, wu) + local(uw, vu, vw)
        + local(wu, vw, vu) + local(vw, uv, uw) + local(wv, uw, uv);
      // i --> j, i --> k, j --> k
      counts(1, 1, 0) += local(uv, uw, vw) + local(vu, vw, uw) + local(uw, uv, wv)
        + local(wu, wv, uv) + local(vw, vu, wu) + local(wv, wu, vu); 
      // i --> j, i --> k, k --> j  int uv = 0, vu = 1, uw = 2, wu = 3, vw = 4, wv = 5;
      counts(1, 1, 1) += local(uv, uw, wv) + local(vu, vw, wu) + local(uw, uv, vw)
        + local(wu, wv, vu) + local(vw, vu, uw) + local(wv, wu, uv);
    }
  }
}

void TempMotifCounter::AddTriadEdgeData(TVec<TriadEdgeData>& events,
                                        TVec<TIntPair>& ts_indices,
                                        int& index, int u, int v, int nbr,
                                        int key1, int key2) {
  if (HasEdges(u, v)) {
    const TIntV& timestamps = temporal_data_[u].GetDat(v);
    for (int i = 0; i < timestamps.Len(); i++) {
      ts_indices.Add(TIntPair(timestamps[i], index));
      events.Add(TriadEdgeData(nbr, key1, key2));
      ++index;
    }
  }
}

//void TempMotifCounter::Count3TEdgeTriads(double delta, Counter3D& counts, TVec<TVec<TStr>>* motifRecordFile) {
//  counts = Counter3D(2, 2, 2);
//
//  // Get the counts on each undirected edge
//  TVec< THash<TInt, TInt> > edge_counts(static_graph_->GetMxNId());
//  TVec< THash<TInt, TIntV> > assignments(static_graph_->GetMxNId());
//  for (TNGraph::TEdgeI it = static_graph_->BegEI();
//       it < static_graph_->EndEI(); it++) {
//    int src = it.GetSrcNId();
//    int dst = it.GetDstNId();
//    int min_node = MIN(src, dst);
//    int max_node = MAX(src, dst);
//    edge_counts[min_node](max_node) += temporal_data_[src](dst).Len();
//    assignments[min_node](max_node) = TIntV();
//  }
//  
//  // Assign triangles to the edge with the most events
//  TIntV Us, Vs, Ws;
//  GetAllStaticTriangles(Us, Vs, Ws);
//  #pragma omp parallel for schedule(dynamic)
//  for (int i = 0; i < Us.Len(); i++) {
//    int u = Us[i];
//    int v = Vs[i];
//    int w = Ws[i];
//    int counts_uv = edge_counts[MIN(u, v)].GetDat(MAX(u, v));
//    int counts_uw = edge_counts[MIN(u, w)].GetDat(MAX(u, w));
//    int counts_vw = edge_counts[MIN(v, w)].GetDat(MAX(v, w));
//    if        (counts_uv >= MAX(counts_uw, counts_vw)) {
//      #pragma omp critical
//      {
//        TIntV& assignment = assignments[MIN(u, v)].GetDat(MAX(u, v));
//        assignment.Add(w);
//      }
//    } else if (counts_uw >= MAX(counts_uv, counts_vw)) {
//      #pragma omp critical
//      {
//        TIntV& assignment = assignments[MIN(u, w)].GetDat(MAX(u, w));
//        assignment.Add(v);      
//      }
//    } else if (counts_vw >= MAX(counts_uv, counts_uw)) {
//      #pragma omp critical
//      {
//        TIntV& assignment = assignments[MIN(v, w)].GetDat(MAX(v, w));
//        assignment.Add(u);              
//      }
//    }
//  }
//
//  TVec<TIntPair> all_edges;
//  TIntV all_nodes;
//  GetAllNodes(all_nodes);  
//  for (int node_id = 0; node_id < all_nodes.Len(); node_id++) {
//    int u = all_nodes[node_id];
//    TIntV nbrs;
//    GetAllNeighbors(u, nbrs);
//    for (int nbr_id = 0; nbr_id < nbrs.Len(); nbr_id++) {
//      int v = nbrs[nbr_id];
//      if (assignments[u].IsKey(v) && assignments[u].GetDat(v).Len() > 0) {
//        all_edges.Add(TIntPair(u, v));
//      }
//    }
//  }
//
//  // Count triangles on edges with the assigned neighbors
//  #pragma omp parallel for schedule(dynamic)
//  for (int edge_id = 0; edge_id < all_edges.Len(); edge_id++) {
//    TIntPair edge = all_edges[edge_id];
//    int u = edge.Key;
//    int v = edge.Dat;
//    // Continue if no assignment
//    if (!assignments[u].IsKey(v)) { continue; }
//    TIntV& uv_assignment = assignments[u].GetDat(v);
//    // Continue if no data
//    if (uv_assignment.Len() == 0) { continue; }
//    // Get all events on (u, v)
//    TVec<TriadEdgeData> events;
//    TVec<TIntPair> ts_indices;
//    int index = 0;
//    int nbr_index = 0;
//    // Assign indices from 0, 1, ..., num_nbrs + 2
//    AddTriadEdgeData(events, ts_indices, index, u, v, nbr_index, 0, 1);
//    nbr_index++;
//    AddTriadEdgeData(events, ts_indices, index, v, u, nbr_index, 0, 0);
//    nbr_index++;
//    // Get all events on triangles assigned to (u, v)
//    for (int w_id = 0; w_id < uv_assignment.Len(); w_id++) {
//      int w = uv_assignment[w_id];
//      AddTriadEdgeData(events, ts_indices, index, w, u, nbr_index, 0, 0);
//      AddTriadEdgeData(events, ts_indices, index, w, v, nbr_index, 0, 1);
//      AddTriadEdgeData(events, ts_indices, index, u, w, nbr_index, 1, 0);
//      AddTriadEdgeData(events, ts_indices, index, v, w, nbr_index, 1, 1);
//      nbr_index++;      
//    }
//    // Put events in sorted order
//    ts_indices.Sort();
//    TIntV timestamps(ts_indices.Len());
//    TVec<TriadEdgeData> sorted_events(ts_indices.Len());
//    TVec<TIntPair> uv(ts_indices.Len());
//    for (int i = 0; i < ts_indices.Len(); i++) {
//      timestamps[i] = ts_indices[i].Key;
//      sorted_events[i] = events[ts_indices[i].Dat];
//
//    }
//    
//    // Get the counts and update the counter
//    ThreeTEdgeTriadCounter tetc(nbr_index, 0, 1);
//    tetc.Count(sorted_events, timestamps, delta, motifRecordFile);
//    #pragma omp critical
//    {
//      for (int dir1 = 0; dir1 < 2; dir1++) {
//        for (int dir2 = 0; dir2 < 2; dir2++) {
//          for (int dir3 = 0; dir3 < 2; dir3++) {        
//            counts(dir1, dir2, dir3) += tetc.Counts(dir1, dir2, dir3);
//          }
//        }
//      }
//    }
//  }
//}

///////////////////////////////////////////////////////////////////////////////
// Generic three temporal edge motif counter
void ThreeTEdgeMotifCounter::Count(const TIntV& event_string, const TIntV& timestamps, const TVec<TIntPair>& uv,
                                   double delta, Counter3D& counts, TVec<TVec<TStr>>* motifRecordFile) {
  // Initialize everything to empty
  counts1_ = Counter1D(size_);
  counts2_ = Counter2D(size_, size_);
  counts3_ = Counter3D(size_, size_, size_);
  //MotifR=TQuad<TInt, TInt, TIntPair, TInt> TRecord=<lable, event, TIntPair uv, TInt timestamps>
  TVec<TRecord> MotifR;

  if (event_string.Len() != timestamps.Len()) {
    TExcept::Throw("Number of events must equal number of timestamps");
  }
  int start = 0;
  for (int end = 0; end < event_string.Len(); end++) {
      while (double(timestamps[start]) + delta < double(timestamps[end])) {
          PrintMotif(&MotifR, motifRecordFile);
          DecrementCounts(event_string[start], uv[start], &MotifR, timestamps[start]);
          start++;
      }
      IncrementCounts(event_string[end], uv[end], &MotifR, timestamps[end], motifRecordFile);
  }
  
  
  counts = counts3_;
 }

//void ThreeTEdgeMotifCounter::EraseSame(TVec<TRecord>* MotifR) {
//    set<int> s(*MotifR.begin(), *MotifR.end());
//    vec.assign(s.begin(), s.end());
//}

void ThreeTEdgeMotifCounter::PrintMotif(TVec<TRecord> *MotifR, TVec<TVec<TStr>>* motifRecordFile) {
    TVec<TRecord> MotifR1;
    TVec<TRecord> MotifR2;
    TVec<TRecord> MotifR3;
    TStr MotifRStr;
    //printf("\nMotifR_len:%d", (*MotifR).Len());
    //Classify MotifR into MotifR1,MotifR2 and MotifR3 according to the lable(0,1,2)
    //TRecord <TIntQu, TIntPair, TInt>:<(lable,i,j,event), uv, timestamps>
    //counts1_:TIntQu(0, 100, 100, event)
    //counts2_:TIntQu(1, i, 100, event)
    //counts3_:TIntQu(2, i, j, event)
    //counts2_(i, event) += counts1_(i)
    //counts3_(i, j, event) += counts2_(i, j)

    for (int i = 0; i < (*MotifR).Len(); i++) {
        if ((*MotifR)[i].Val1.Val1 == 0) {
            MotifR1.Add((*MotifR)[i]);
        }
        else if ((*MotifR)[i].Val1.Val1 == 1) {
            MotifR2.Add((*MotifR)[i]);
        }
        else {
            MotifR3.Add((*MotifR)[i]);
        }
    }
    /*for (int i = 0; i < (*MotifR).Len(); i++) {
        printf("\nMotifR:%d-%d-%d-%d-%d", (*MotifR)[i].Val1, (*MotifR)[i].Val2, (*MotifR)[i].Val3.Key, (*MotifR)[i].Val3.Dat, (*MotifR)[i].Val4);
    }*/

    //*MotifR=TRecord <TIntQu, TIntPair, TInt>:<(lable,i,j,event), uv, timestamps>
    
    for (int i = 0; i < MotifR1.Len(); i++) {
        for (int j = 0; j < MotifR2.Len(); j++) {
            for (int l = 0; l < MotifR3.Len(); l++) {
                //2 Node Motif 
                if ((MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (5 - 1) * 6 + 1 - 1;
                    MotifRStr = "M_{5,1}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                {
                    int index = 36 + 2 * 2 * 1 + 2 * 0 + 1 * 1;
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr());}
                }
                else if (MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                {
                    int index = 36 + 2 * 2 * 0 + 2 * 1 + 1 * 0;
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }

                else if ((MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (5 - 1) * 6 + 2 - 1;
                    MotifRStr = "M_{5,2}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                {
                    int index = 36 + 2 * 2 * 1 + 2 * 0 + 1 * 0;
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                {
                    int index = 36 + 2 * 2 * 0 + 2 * 1 + 1 * 1;
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (6 - 1) * 6 + 1 - 1;
                    MotifRStr = "M_{6,1}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                {
                    int index = 36 + 2 * 2 * 0 + 2 * 0 + 1 * 0;
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                {
                    int index = 36 + 2 * 2 * 1 + 2 * 1 + 1 * 1;
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (6 - 1) * 6 + 2 - 1;
                    MotifRStr = "M_{6,2}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                {
                    int index = 36 + 2 * 2 * 0 + 2 * 0 + 1 * 1;
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                {
                    int index = 36 + 2 * 2 * 1 + 2 * 1 + 1 * 0;
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }

                //2 Node Motif 
                /*if ((MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (5 - 1) * 6 + 1 - 1;
                    MotifRStr = "M_{5,1}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (5 - 1) * 6 + 2 - 1;
                    MotifRStr = "M_{5,2}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (6 - 1) * 6 + 1 - 1;
                    MotifRStr = "M_{6,1}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (6 - 1) * 6 + 2 - 1;
                    MotifRStr = "M_{6,2}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }*/

                 //3 Node Star Motif
                //value restrict
                //time restrict
                //index value restrict

                else if ((MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 3 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 3 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 3 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (1 - 1) * 6 + 1 - 1;
                    MotifRStr = "M_{1,1}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 3 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 3 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 2 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (1 - 1) * 6 + 2 - 1;
                    MotifRStr = "M_{1,2}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 3 && MotifR3[l].Val1.Val4 == 2 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 3 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (1 - 1) * 6 + 5 - 1;
                    MotifRStr = "M_{1,5}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 3 && MotifR3[l].Val1.Val4 == 3 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 3 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (1 - 1) * 6 + 6 - 1;
                    MotifRStr = "M_{1,6}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 2 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 3 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 3 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (2 - 1) * 6 + 1 - 1;
                    MotifRStr = "M_{2,1}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 2 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 3 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 2 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (2 - 1) * 6 + 2 - 1;
                    MotifRStr = "M_{2,2}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 2 && MotifR3[l].Val1.Val4 == 2 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 3 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (2 - 1) * 6 + 5 - 1;
                    MotifRStr = "M_{2,5}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 2 && MotifR3[l].Val1.Val4 == 3 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 3 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (2 - 1) * 6 + 6 - 1;
                    MotifRStr = "M_{2,6}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 3 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 2 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 2 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (3 - 1) * 6 + 1 - 1;
                    MotifRStr = "M_{3,1}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 3 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 2 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 3 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (3 - 1) * 6 + 2 - 1;
                    MotifRStr = "M_{3,2}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 3 && MotifR3[l].Val1.Val4 == 2 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 2 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (3 - 1) * 6 + 3 - 1;
                    MotifRStr = "M_{3,3}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 3 && MotifR3[l].Val1.Val4 == 3 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 2 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (3 - 1) * 6 + 4 - 1;
                    MotifRStr = "M_{3,4}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 2 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 2 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 2 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (4 - 1) * 6 + 1 - 1;
                    MotifRStr = "M_{4,1}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 2 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 2 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 3 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (4 - 1) * 6 + 2 - 1;
                    MotifRStr = "M_{4,2}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 2 && MotifR3[l].Val1.Val4 == 2 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 2 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (4 - 1) * 6 + 3 - 1;
                    MotifRStr = "M_{4,3}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 2 && MotifR3[l].Val1.Val4 == 3 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 2 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (4 - 1) * 6 + 4 - 1;
                    MotifRStr = "M_{4,4}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 2 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 2 && MotifR3[l].Val1.Val3 == 3 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (5 - 1) * 6 + 3 - 1;
                    MotifRStr = "M_{5,3}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 3 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 2 && MotifR3[l].Val1.Val3 == 3 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (5 - 1) * 6 + 4 - 1;
                    MotifRStr = "M_{5,4}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 2 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 3 && MotifR3[l].Val1.Val3 == 2 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (5 - 1) * 6 + 5 - 1;
                    MotifRStr = "M_{5,5}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 3 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 3 && MotifR3[l].Val1.Val3 == 2 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (5 - 1) * 6 + 6 - 1;
                    MotifRStr = "M_{5,6}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 2 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 2 && MotifR3[l].Val1.Val3 == 2 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (6 - 1) * 6 + 3 - 1;
                    MotifRStr = "M_{6,3}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 3 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 2 && MotifR3[l].Val1.Val3 == 2 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (6 - 1) * 6 + 4 - 1;
                    MotifRStr = "M_{6,4}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 2 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 3 && MotifR3[l].Val1.Val3 == 3 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (6 - 1) * 6 + 5 - 1;
                    MotifRStr = "M_{6,5}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if ((MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 3 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                    ||
                    (MotifR3[l].Val1.Val2 == 3 && MotifR3[l].Val1.Val3 == 3 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                {
                    int index = (6 - 1) * 6 + 6 - 1;
                    MotifRStr = "M_{6,6}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }


                 //3 Node Triad Motif
                 else if ((MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 5 && MotifR3[l].Val1.Val4 == 2 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 3 && MotifR3[l].Val1.Val4 == 4 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 2 && MotifR3[l].Val1.Val3 == 4 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 3 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 5 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 4 && MotifR3[l].Val1.Val3 == 2 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 5 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 3 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                 {
                     int index = (1 - 1) * 6 + 3 - 1;
                     MotifRStr = "M_{1,3}:\n";
                     if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                     MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                     if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                 }
                 else if ((MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 5 && MotifR3[l].Val1.Val4 == 3 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 3 && MotifR3[l].Val1.Val4 == 5 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 2 && MotifR3[l].Val1.Val3 == 4 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 3 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 4 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 4 && MotifR3[l].Val1.Val3 == 2 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 5 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 2 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                 {
                     int index = (1 - 1) * 6 + 4 - 1;
                     MotifRStr = "M_{1,4}:\n";
                     if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                     MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                     if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                 }
                 else if ((MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 4 && MotifR3[l].Val1.Val4 == 2 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 2 && MotifR3[l].Val1.Val4 == 4 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 2 && MotifR3[l].Val1.Val3 == 5 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 3 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 5 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 4 && MotifR3[l].Val1.Val3 == 3 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 5 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 3 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                 {
                     int index = (2 - 1) * 6 + 3 - 1;
                     MotifRStr = "M_{2,3}:\n";
                     if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                     MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                     if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                 }
                 else if ((MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 4 && MotifR3[l].Val1.Val4 == 3 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 2 && MotifR3[l].Val1.Val4 == 5 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 2 && MotifR3[l].Val1.Val3 == 5 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 3 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 4 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 4 && MotifR3[l].Val1.Val3 == 3 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 5 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 2 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                 {
                     int index = (2 - 1) * 6 + 4 - 1;
                     MotifRStr = "M_{2,4}:\n";
                     if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                     MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                     if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                 }
                else if ((MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 5 && MotifR3[l].Val1.Val4 == 2 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 3 && MotifR3[l].Val1.Val4 == 4 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 3 && MotifR3[l].Val1.Val3 == 4 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 2 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 5 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 5 && MotifR3[l].Val1.Val3 == 2 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 4 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 3 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                 {
                     int index = (3 - 1) * 6 + 5 - 1;
                     MotifRStr = "M_{3,5}:\n";
                     if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                     MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                     if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                 }
                else if ((MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 5 && MotifR3[l].Val1.Val4 == 3 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 3 && MotifR3[l].Val1.Val4 == 5 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 3 && MotifR3[l].Val1.Val3 == 4 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 2 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 4 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 5 && MotifR3[l].Val1.Val3 == 2 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 4 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 2 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                 {
                     int index = (3 - 1) * 6 + 6 - 1;
                     MotifRStr = "M_{3,6}:\n";
                     if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                     MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                     if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                 }
                 else if ((MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 4 && MotifR3[l].Val1.Val4 == 2 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 2 && MotifR3[l].Val1.Val4 == 4 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 3 && MotifR3[l].Val1.Val3 == 5 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 2 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 5 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 5 && MotifR3[l].Val1.Val3 == 3 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 4 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 3 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                 {
                     int index = (4 - 1) * 6 + 5 - 1;
                     MotifRStr = "M_{4,5}:\n";
                     if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                     MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                     if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                 }
                 else if ((MotifR3[l].Val1.Val2 == 1 && MotifR3[l].Val1.Val3 == 4 && MotifR3[l].Val1.Val4 == 3 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 0 && MotifR3[l].Val1.Val3 == 2 && MotifR3[l].Val1.Val4 == 5 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 3 && MotifR3[l].Val1.Val3 == 5 && MotifR3[l].Val1.Val4 == 1 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 2 && MotifR3[l].Val1.Val3 == 0 && MotifR3[l].Val1.Val4 == 4 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 5 && MotifR3[l].Val1.Val3 == 3 && MotifR3[l].Val1.Val4 == 0 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3)
                     ||
                     (MotifR3[l].Val1.Val2 == 4 && MotifR3[l].Val1.Val3 == 1 && MotifR3[l].Val1.Val4 == 2 && MotifR1[i].Val3 < MotifR2[j].Val3 && MotifR2[j].Val3 < MotifR3[l].Val3 && MotifR1[i].Val1.Val4 == MotifR2[j].Val1.Val2 && MotifR1[i].Val1.Val4 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val2 == MotifR3[l].Val1.Val2 && MotifR2[j].Val1.Val4 == MotifR3[l].Val1.Val3))
                 {
                     int index = (4 - 1) * 6 + 6 - 1;
                     MotifRStr = "M_{4,6}:\n";
                     if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                     MotifRStr = "(" + MotifR1[i].Val2.Key.GetStr() + "," + MotifR1[i].Val2.Dat.GetStr() + "," + MotifR1[i].Val3.GetStr() + ")+(" + MotifR2[j].Val2.Key.GetStr() + "," + MotifR2[j].Val2.Dat.GetStr() + "," + MotifR2[j].Val3.GetStr() + ")+(" + MotifR3[l].Val2.Key.GetStr() + "," + MotifR3[l].Val2.Dat.GetStr() + "," + MotifR3[l].Val3.GetStr() + ")" + "\n";
                     if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                 }
            }
        }
        
    } 
}

void ThreeTEdgeMotifCounter::IncrementCounts(int event, TIntPair uv, TVec<TRecord> *MotifR, TInt timestamps, TVec<TVec<TStr>>* motifRecordFile) {
  //printf("\ninitial MotifR:%d", MotifR.Len());
  for (int i = 0; i < size_; i++) {
    for (int j = 0; j < size_; j++) {
      int k = counts3_(i, j, event);
      counts3_(i, j, event) += counts2_(i, j);
      if (k < counts3_(i, j, event)&& !(*MotifR).IsIn(TRecord(TIntQu(2, i, j, event), uv, timestamps))) {
          //printf("\nNew MotifR:%d", counts2_(i, j));
          (*MotifR).Add(TRecord(TIntQu(2, i, j, event), uv, timestamps));
          //PrintMotif(MotifR, motifRecordFile);
      }
    }
  }
  // 100 can be replaced by any number bigger than 5, in order to be different from _size value
  for (int i = 0; i < size_; i++) {
      int k = counts2_(i, event);
      counts2_(i, event) += counts1_(i);
      if (k < counts2_(i, event) && !(*MotifR).IsIn(TRecord(TIntQu(1, i, 100, event), uv, timestamps))) {
          (*MotifR).Add(TRecord(TIntQu(1, i, 100, event), uv, timestamps));
      }
  }
  counts1_(event) += 1;
  if (!(*MotifR).IsIn(TRecord(TIntQu(0, 100, 100, event), uv, timestamps))) {
      (*MotifR).Add(TRecord(TIntQu(0, 100, 100, event), uv, timestamps));
  }
}

void ThreeTEdgeMotifCounter::DecrementCounts(int event, TIntPair uv, TVec<TRecord> *MotifR, TInt timestamps) {
    //printf("\ninitial MotifR:%d", MotifR.Len());
    counts1_(event)--;
    (*MotifR).DelAll(TRecord(TIntQu(0, 100, 100, event), uv, timestamps));
    for (int i = 0; i < size_; i++) {
        int k = counts2_(i, event);
        counts2_(event, i) -= counts1_(i);
        if (k > counts2_(i, event)) {
            (*MotifR).DelAll(TRecord(TIntQu(1, i, 100, event), uv, timestamps));
        } 
    }
  //printf("\nDecre MotifR:%d", (*MotifR).Len());
}


//void ThreeTEdgeMotifCounter::Count(const TIntV& event_string, const TIntV& timestamps,
//    double delta, Counter3D& counts) {
//    // Initialize everything to empty
//    counts1_ = Counter1D(size_);
//    counts2_ = Counter2D(size_, size_);
//    counts3_ = Counter3D(size_, size_, size_);
//    //cout << "size:" << size_ << endl;
//    if (event_string.Len() != timestamps.Len()) {
//        TExcept::Throw("Number of events must equal number of timestamps");
//    }
//    int start = 0;
//    //printf("\nevent_string:%d", event_string.Len());
//    //printf("\ntimestamps:%d", timestamps.Len());
//    for (int end = 0; end < event_string.Len(); end++) {
//        while (double(timestamps[start]) + delta < double(timestamps[end])) {
//            DecrementCounts(event_string[start]);
//            start++;
//        }
//        IncrementCounts(event_string[end]);
//    }
//    counts = counts3_;
//}
//
//void ThreeTEdgeMotifCounter::IncrementCounts(int event) {
//    for (int i = 0; i < size_; i++) {
//        for (int j = 0; j < size_; j++) {
//            counts3_(i, j, event) += counts2_(i, j);
//        }
//    }
//    for (int i = 0; i < size_; i++) { counts2_(i, event) += counts1_(i); }
//    counts1_(event) += 1;
//}
//
//void ThreeTEdgeMotifCounter::DecrementCounts(int event) {
//    counts1_(event)--;
//    for (int i = 0; i < size_; i++) { counts2_(event, i) -= counts1_(i); }
//}

void ThreeTEdgeStarCounter::PrintStarMotif(TVec<TRecordStar>* MotifRStar, TVec<TVec<TStr>>* motifRecordFile) {
    TVec<TRecordStar> MotifRStar1;
    TVec<TRecordStar> MotifRStar2;
    TVec<TRecordStar> MotifRStar3;
    TStr MotifRStr;
    //*MotifRStar=TRecordStar <TIntQu, TIntPair, TInt>:<(level,i,j,event), uv, (lable,timestamps)>
    for (int i = 0; i < (*MotifRStar).Len(); i++) {
        if ((*MotifRStar)[i].Val1.Val1 == 0) {
            MotifRStar1.Add((*MotifRStar)[i]);
        }
        else if ((*MotifRStar)[i].Val1.Val1 == 1) {
            MotifRStar2.Add((*MotifRStar)[i]);
        }
        else {
            MotifRStar3.Add((*MotifRStar)[i]);
        }
    }

    for (int i = 0; i < MotifRStar1.Len(); i++) {
        for (int j = 0; j < MotifRStar2.Len(); j++) {
            for (int l = 0; l < MotifRStar3.Len(); l++) {
                //pre=0  mid=1  pos=2
                //count=3 
                //sum=2
                //node=1
                //time sequence

                //mid
                if (MotifRStar1[i].Val3.Key == 2 && MotifRStar2[j].Val3.Key == 1 && MotifRStar3[l].Val3.Key == 1 && 
                    MotifRStar3[l].Val1.Val2 == 1 && MotifRStar3[l].Val1.Val3 == 1 && MotifRStar3[l].Val1.Val4 == 1 && 
                    MotifRStar2[j].Val1.Val2 == 1 && MotifRStar2[j].Val1.Val3 == 1 && 
                    MotifRStar1[i].Val1.Val2 == 1 &&
                    MotifRStar1[i].Val3.Dat < MotifRStar2[j].Val3.Dat && MotifRStar2[j].Val3.Dat < MotifRStar3[l].Val3.Dat)
                {
                    int index = (1 - 1) * 6 + 1 - 1;
                    MotifRStr = "M_{1,1}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifRStar1[i].Val2.Key.GetStr() + "," + MotifRStar1[i].Val2.Dat.GetStr() + "," + MotifRStar1[i].Val3.Dat.GetStr() + ")+(" + MotifRStar2[j].Val2.Key.GetStr() + "," + MotifRStar2[j].Val2.Dat.GetStr() + "," + MotifRStar2[j].Val3.Dat.GetStr() + ")+(" + MotifRStar3[l].Val2.Key.GetStr() + "," + MotifRStar3[l].Val2.Dat.GetStr() + "," + MotifRStar3[l].Val3.Dat.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifRStar1[i].Val3.Key == 2 && MotifRStar2[j].Val3.Key == 1 && MotifRStar3[l].Val3.Key == 1 &&
                    MotifRStar3[l].Val1.Val2 == 1 && MotifRStar3[l].Val1.Val3 == 1 && MotifRStar3[l].Val1.Val4 == 0 &&
                    MotifRStar2[j].Val1.Val2 == 0 && MotifRStar2[j].Val1.Val3 == 1 &&
                    MotifRStar1[i].Val1.Val2 == 0 &&
                    MotifRStar1[i].Val3.Dat < MotifRStar2[j].Val3.Dat && MotifRStar2[j].Val3.Dat < MotifRStar3[l].Val3.Dat)
                {
                    int index = (1 - 1) * 6 + 2 - 1;
                    MotifRStr = "M_{1,2}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifRStar1[i].Val2.Key.GetStr() + "," + MotifRStar1[i].Val2.Dat.GetStr() + "," + MotifRStar1[i].Val3.Dat.GetStr() + ")+(" + MotifRStar2[j].Val2.Key.GetStr() + "," + MotifRStar2[j].Val2.Dat.GetStr() + "," + MotifRStar2[j].Val3.Dat.GetStr() + ")+(" + MotifRStar3[l].Val2.Key.GetStr() + "," + MotifRStar3[l].Val2.Dat.GetStr() + "," + MotifRStar3[l].Val3.Dat.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifRStar1[i].Val3.Key == 2 && MotifRStar2[j].Val3.Key == 1 && MotifRStar3[l].Val3.Key == 1 &&
                    MotifRStar3[l].Val1.Val2 == 1 && MotifRStar3[l].Val1.Val3 == 0 && MotifRStar3[l].Val1.Val4 == 1 &&
                    MotifRStar2[j].Val1.Val2 == 1 && MotifRStar2[j].Val1.Val3 == 1 &&
                    MotifRStar1[i].Val1.Val2 == 1 &&
                    MotifRStar1[i].Val3.Dat < MotifRStar2[j].Val3.Dat && MotifRStar2[j].Val3.Dat < MotifRStar3[l].Val3.Dat)
                {
                    int index = (2 - 1) * 6 + 1 - 1;
                    MotifRStr = "M_{2,1}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifRStar1[i].Val2.Key.GetStr() + "," + MotifRStar1[i].Val2.Dat.GetStr() + "," + MotifRStar1[i].Val3.Dat.GetStr() + ")+(" + MotifRStar2[j].Val2.Key.GetStr() + "," + MotifRStar2[j].Val2.Dat.GetStr() + "," + MotifRStar2[j].Val3.Dat.GetStr() + ")+(" + MotifRStar3[l].Val2.Key.GetStr() + "," + MotifRStar3[l].Val2.Dat.GetStr() + "," + MotifRStar3[l].Val3.Dat.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifRStar1[i].Val3.Key == 2 && MotifRStar2[j].Val3.Key == 1 && MotifRStar3[l].Val3.Key == 1 &&
                    MotifRStar3[l].Val1.Val2 == 1 && MotifRStar3[l].Val1.Val3 == 0 && MotifRStar3[l].Val1.Val4 == 0 &&
                    MotifRStar2[j].Val1.Val2 == 0 && MotifRStar2[j].Val1.Val3 == 1 &&
                    MotifRStar1[i].Val1.Val2 == 0 &&
                    MotifRStar1[i].Val3.Dat < MotifRStar2[j].Val3.Dat && MotifRStar2[j].Val3.Dat < MotifRStar3[l].Val3.Dat)
                {
                    int index = (2 - 1) * 6 + 2 - 1;
                    MotifRStr = "M_{2,2}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifRStar1[i].Val2.Key.GetStr() + "," + MotifRStar1[i].Val2.Dat.GetStr() + "," + MotifRStar1[i].Val3.Dat.GetStr() + ")+(" + MotifRStar2[j].Val2.Key.GetStr() + "," + MotifRStar2[j].Val2.Dat.GetStr() + "," + MotifRStar2[j].Val3.Dat.GetStr() + ")+(" + MotifRStar3[l].Val2.Key.GetStr() + "," + MotifRStar3[l].Val2.Dat.GetStr() + "," + MotifRStar3[l].Val3.Dat.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifRStar1[i].Val3.Key == 2 && MotifRStar2[j].Val3.Key == 1 && MotifRStar3[l].Val3.Key == 1 &&
                    MotifRStar3[l].Val1.Val2 == 0 && MotifRStar3[l].Val1.Val3 == 1 && MotifRStar3[l].Val1.Val4 == 0 &&
                    MotifRStar2[j].Val1.Val2 == 0 && MotifRStar2[j].Val1.Val3 == 0 &&
                    MotifRStar1[i].Val1.Val2 == 0 &&
                    MotifRStar1[i].Val3.Dat < MotifRStar2[j].Val3.Dat && MotifRStar2[j].Val3.Dat < MotifRStar3[l].Val3.Dat)
                {
                    int index = (3 - 1) * 6 + 1 - 1;
                    MotifRStr = "M_{3,1}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifRStar1[i].Val2.Key.GetStr() + "," + MotifRStar1[i].Val2.Dat.GetStr() + "," + MotifRStar1[i].Val3.Dat.GetStr() + ")+(" + MotifRStar2[j].Val2.Key.GetStr() + "," + MotifRStar2[j].Val2.Dat.GetStr() + "," + MotifRStar2[j].Val3.Dat.GetStr() + ")+(" + MotifRStar3[l].Val2.Key.GetStr() + "," + MotifRStar3[l].Val2.Dat.GetStr() + "," + MotifRStar3[l].Val3.Dat.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifRStar1[i].Val3.Key == 2 && MotifRStar2[j].Val3.Key == 1 && MotifRStar3[l].Val3.Key == 1 &&
                    MotifRStar3[l].Val1.Val2 == 0 && MotifRStar3[l].Val1.Val3 == 1 && MotifRStar3[l].Val1.Val4 == 1 &&
                    MotifRStar2[j].Val1.Val2 == 1 && MotifRStar2[j].Val1.Val3 == 0 &&
                    MotifRStar1[i].Val1.Val2 == 1 &&
                    MotifRStar1[i].Val3.Dat < MotifRStar2[j].Val3.Dat && MotifRStar2[j].Val3.Dat < MotifRStar3[l].Val3.Dat)
                {
                    int index = (3 - 1) * 6 + 2 - 1;
                    MotifRStr = "M_{3,2}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifRStar1[i].Val2.Key.GetStr() + "," + MotifRStar1[i].Val2.Dat.GetStr() + "," + MotifRStar1[i].Val3.Dat.GetStr() + ")+(" + MotifRStar2[j].Val2.Key.GetStr() + "," + MotifRStar2[j].Val2.Dat.GetStr() + "," + MotifRStar2[j].Val3.Dat.GetStr() + ")+(" + MotifRStar3[l].Val2.Key.GetStr() + "," + MotifRStar3[l].Val2.Dat.GetStr() + "," + MotifRStar3[l].Val3.Dat.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifRStar1[i].Val3.Key == 2 && MotifRStar2[j].Val3.Key == 1 && MotifRStar3[l].Val3.Key == 1 &&
                    MotifRStar3[l].Val1.Val2 == 0 && MotifRStar3[l].Val1.Val3 == 0 && MotifRStar3[l].Val1.Val4 == 0 &&
                    MotifRStar2[j].Val1.Val2 == 0 && MotifRStar2[j].Val1.Val3 == 0 &&
                    MotifRStar1[i].Val1.Val2 == 0 &&
                    MotifRStar1[i].Val3.Dat < MotifRStar2[j].Val3.Dat && MotifRStar2[j].Val3.Dat < MotifRStar3[l].Val3.Dat)
                {
                    int index = (4 - 1) * 6 + 1 - 1;
                    MotifRStr = "M_{4,1}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifRStar1[i].Val2.Key.GetStr() + "," + MotifRStar1[i].Val2.Dat.GetStr() + "," + MotifRStar1[i].Val3.Dat.GetStr() + ")+(" + MotifRStar2[j].Val2.Key.GetStr() + "," + MotifRStar2[j].Val2.Dat.GetStr() + "," + MotifRStar2[j].Val3.Dat.GetStr() + ")+(" + MotifRStar3[l].Val2.Key.GetStr() + "," + MotifRStar3[l].Val2.Dat.GetStr() + "," + MotifRStar3[l].Val3.Dat.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifRStar1[i].Val3.Key == 2 && MotifRStar2[j].Val3.Key == 1 && MotifRStar3[l].Val3.Key == 1 &&
                    MotifRStar3[l].Val1.Val2 == 0 && MotifRStar3[l].Val1.Val3 == 0 && MotifRStar3[l].Val1.Val4 == 1 &&
                    MotifRStar2[j].Val1.Val2 == 1 && MotifRStar2[j].Val1.Val3 == 0 &&
                    MotifRStar1[i].Val1.Val2 == 1 &&
                    MotifRStar1[i].Val3.Dat < MotifRStar2[j].Val3.Dat && MotifRStar2[j].Val3.Dat < MotifRStar3[l].Val3.Dat)
                {
                    int index = (3 - 1) * 6 + 2 - 1;
                    MotifRStr = "M_{4,2}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifRStar1[i].Val2.Key.GetStr() + "," + MotifRStar1[i].Val2.Dat.GetStr() + "," + MotifRStar1[i].Val3.Dat.GetStr() + ")+(" + MotifRStar2[j].Val2.Key.GetStr() + "," + MotifRStar2[j].Val2.Dat.GetStr() + "," + MotifRStar2[j].Val3.Dat.GetStr() + ")+(" + MotifRStar3[l].Val2.Key.GetStr() + "," + MotifRStar3[l].Val2.Dat.GetStr() + "," + MotifRStar3[l].Val3.Dat.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                //pos
                else if (MotifRStar1[i].Val3.Key == 2 && MotifRStar2[j].Val3.Key == 2 && MotifRStar3[l].Val3.Key == 2 &&
                    MotifRStar3[l].Val1.Val2 == 1 && MotifRStar3[l].Val1.Val3 == 1 && MotifRStar3[l].Val1.Val4 == 0 &&
                    MotifRStar2[j].Val1.Val2 == 1 && MotifRStar2[j].Val1.Val3 == 0 &&
                    MotifRStar1[i].Val1.Val2 == 1 &&
                    MotifRStar1[i].Val3.Dat < MotifRStar2[j].Val3.Dat && MotifRStar2[j].Val3.Dat < MotifRStar3[l].Val3.Dat)
                {
                    int index = (1 - 1) * 6 + 5 - 1;
                    MotifRStr = "M_{1,5}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifRStar1[i].Val2.Key.GetStr() + "," + MotifRStar1[i].Val2.Dat.GetStr() + "," + MotifRStar1[i].Val3.Dat.GetStr() + ")+(" + MotifRStar2[j].Val2.Key.GetStr() + "," + MotifRStar2[j].Val2.Dat.GetStr() + "," + MotifRStar2[j].Val3.Dat.GetStr() + ")+(" + MotifRStar3[l].Val2.Key.GetStr() + "," + MotifRStar3[l].Val2.Dat.GetStr() + "," + MotifRStar3[l].Val3.Dat.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifRStar1[i].Val3.Key == 2 && MotifRStar2[j].Val3.Key == 2 && MotifRStar3[l].Val3.Key == 2 &&
                    MotifRStar3[l].Val1.Val2 == 1 && MotifRStar3[l].Val1.Val3 == 1 && MotifRStar3[l].Val1.Val4 == 1 &&
                    MotifRStar2[j].Val1.Val2 == 1 && MotifRStar2[j].Val1.Val3 == 1 &&
                    MotifRStar1[i].Val1.Val2 == 1 &&
                    MotifRStar1[i].Val3.Dat < MotifRStar2[j].Val3.Dat && MotifRStar2[j].Val3.Dat < MotifRStar3[l].Val3.Dat)
                {
                    int index = (1 - 1) * 6 + 6 - 1;
                    MotifRStr = "M_{1,6}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifRStar1[i].Val2.Key.GetStr() + "," + MotifRStar1[i].Val2.Dat.GetStr() + "," + MotifRStar1[i].Val3.Dat.GetStr() + ")+(" + MotifRStar2[j].Val2.Key.GetStr() + "," + MotifRStar2[j].Val2.Dat.GetStr() + "," + MotifRStar2[j].Val3.Dat.GetStr() + ")+(" + MotifRStar3[l].Val2.Key.GetStr() + "," + MotifRStar3[l].Val2.Dat.GetStr() + "," + MotifRStar3[l].Val3.Dat.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifRStar1[i].Val3.Key == 2 && MotifRStar2[j].Val3.Key == 2 && MotifRStar3[l].Val3.Key == 2 &&
                    MotifRStar3[l].Val1.Val2 == 1 && MotifRStar3[l].Val1.Val3 == 0 && MotifRStar3[l].Val1.Val4 == 0 &&
                    MotifRStar2[j].Val1.Val2 == 0 && MotifRStar2[j].Val1.Val3 == 0 &&
                    MotifRStar1[i].Val1.Val2 == 0 &&
                    MotifRStar1[i].Val3.Dat < MotifRStar2[j].Val3.Dat && MotifRStar2[j].Val3.Dat < MotifRStar3[l].Val3.Dat)
                {
                    int index = (2 - 1) * 6 + 5 - 1;
                    MotifRStr = "M_{2,5}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifRStar1[i].Val2.Key.GetStr() + "," + MotifRStar1[i].Val2.Dat.GetStr() + "," + MotifRStar1[i].Val3.Dat.GetStr() + ")+(" + MotifRStar2[j].Val2.Key.GetStr() + "," + MotifRStar2[j].Val2.Dat.GetStr() + "," + MotifRStar2[j].Val3.Dat.GetStr() + ")+(" + MotifRStar3[l].Val2.Key.GetStr() + "," + MotifRStar3[l].Val2.Dat.GetStr() + "," + MotifRStar3[l].Val3.Dat.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifRStar1[i].Val3.Key == 2 && MotifRStar2[j].Val3.Key == 2 && MotifRStar3[l].Val3.Key == 2 &&
                    MotifRStar3[l].Val1.Val2 == 1 && MotifRStar3[l].Val1.Val3 == 0 && MotifRStar3[l].Val1.Val4 == 1 &&
                    MotifRStar2[j].Val1.Val2 == 0 && MotifRStar2[j].Val1.Val3 == 1 &&
                    MotifRStar1[i].Val1.Val2 == 1 &&
                    MotifRStar1[i].Val3.Dat < MotifRStar2[j].Val3.Dat && MotifRStar2[j].Val3.Dat < MotifRStar3[l].Val3.Dat)
                {
                    int index = (2 - 1) * 6 + 6 - 1;
                    MotifRStr = "M_{2,6}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifRStar1[i].Val2.Key.GetStr() + "," + MotifRStar1[i].Val2.Dat.GetStr() + "," + MotifRStar1[i].Val3.Dat.GetStr() + ")+(" + MotifRStar2[j].Val2.Key.GetStr() + "," + MotifRStar2[j].Val2.Dat.GetStr() + "," + MotifRStar2[j].Val3.Dat.GetStr() + ")+(" + MotifRStar3[l].Val2.Key.GetStr() + "," + MotifRStar3[l].Val2.Dat.GetStr() + "," + MotifRStar3[l].Val3.Dat.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifRStar1[i].Val3.Key == 2 && MotifRStar2[j].Val3.Key == 2 && MotifRStar3[l].Val3.Key == 2 &&
                    MotifRStar3[l].Val1.Val2 == 0 && MotifRStar3[l].Val1.Val3 == 1 && MotifRStar3[l].Val1.Val4 == 0 &&
                    MotifRStar2[j].Val1.Val2 == 1 && MotifRStar2[j].Val1.Val3 == 0 &&
                    MotifRStar1[i].Val1.Val2 == 1 &&
                    MotifRStar1[i].Val3.Dat < MotifRStar2[j].Val3.Dat && MotifRStar2[j].Val3.Dat < MotifRStar3[l].Val3.Dat)
                {
                    int index = (3 - 1) * 6 + 3 - 1;
                    MotifRStr = "M_{3,3}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifRStar1[i].Val2.Key.GetStr() + "," + MotifRStar1[i].Val2.Dat.GetStr() + "," + MotifRStar1[i].Val3.Dat.GetStr() + ")+(" + MotifRStar2[j].Val2.Key.GetStr() + "," + MotifRStar2[j].Val2.Dat.GetStr() + "," + MotifRStar2[j].Val3.Dat.GetStr() + ")+(" + MotifRStar3[l].Val2.Key.GetStr() + "," + MotifRStar3[l].Val2.Dat.GetStr() + "," + MotifRStar3[l].Val3.Dat.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifRStar1[i].Val3.Key == 2 && MotifRStar2[j].Val3.Key == 2 && MotifRStar3[l].Val3.Key == 2 &&
                    MotifRStar3[l].Val1.Val2 == 0 && MotifRStar3[l].Val1.Val3 == 1 && MotifRStar3[l].Val1.Val4 == 1 &&
                    MotifRStar2[j].Val1.Val2 == 1 && MotifRStar2[j].Val1.Val3 == 1 &&
                    MotifRStar1[i].Val1.Val2 == 1 &&
                    MotifRStar1[i].Val3.Dat < MotifRStar2[j].Val3.Dat && MotifRStar2[j].Val3.Dat < MotifRStar3[l].Val3.Dat)
                {
                    int index = (3 - 1) * 6 + 4 - 1;
                    MotifRStr = "M_{3,4}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifRStar1[i].Val2.Key.GetStr() + "," + MotifRStar1[i].Val2.Dat.GetStr() + "," + MotifRStar1[i].Val3.Dat.GetStr() + ")+(" + MotifRStar2[j].Val2.Key.GetStr() + "," + MotifRStar2[j].Val2.Dat.GetStr() + "," + MotifRStar2[j].Val3.Dat.GetStr() + ")+(" + MotifRStar3[l].Val2.Key.GetStr() + "," + MotifRStar3[l].Val2.Dat.GetStr() + "," + MotifRStar3[l].Val3.Dat.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifRStar1[i].Val3.Key == 2 && MotifRStar2[j].Val3.Key == 2 && MotifRStar3[l].Val3.Key == 2 &&
                    MotifRStar3[l].Val1.Val2 == 0 && MotifRStar3[l].Val1.Val3 == 0 && MotifRStar3[l].Val1.Val4 == 0 &&
                    MotifRStar2[j].Val1.Val2 == 0 && MotifRStar2[j].Val1.Val3 == 0 &&
                    MotifRStar1[i].Val1.Val2 == 0 &&
                    MotifRStar1[i].Val3.Dat < MotifRStar2[j].Val3.Dat && MotifRStar2[j].Val3.Dat < MotifRStar3[l].Val3.Dat)
                {
                    int index = (4 - 1) * 6 + 3 - 1;
                    MotifRStr = "M_{4,3}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifRStar1[i].Val2.Key.GetStr() + "," + MotifRStar1[i].Val2.Dat.GetStr() + "," + MotifRStar1[i].Val3.Dat.GetStr() + ")+(" + MotifRStar2[j].Val2.Key.GetStr() + "," + MotifRStar2[j].Val2.Dat.GetStr() + "," + MotifRStar2[j].Val3.Dat.GetStr() + ")+(" + MotifRStar3[l].Val2.Key.GetStr() + "," + MotifRStar3[l].Val2.Dat.GetStr() + "," + MotifRStar3[l].Val3.Dat.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifRStar1[i].Val3.Key == 2 && MotifRStar2[j].Val3.Key == 2 && MotifRStar3[l].Val3.Key == 2 &&
                    MotifRStar3[l].Val1.Val2 == 0 && MotifRStar3[l].Val1.Val3 == 0 && MotifRStar3[l].Val1.Val4 == 1 &&
                    MotifRStar2[j].Val1.Val2 == 0 && MotifRStar2[j].Val1.Val3 == 1 &&
                    MotifRStar1[i].Val1.Val2 == 0 &&
                    MotifRStar1[i].Val3.Dat < MotifRStar2[j].Val3.Dat && MotifRStar2[j].Val3.Dat < MotifRStar3[l].Val3.Dat)
                {
                    int index = (4 - 1) * 6 + 4 - 1;
                    MotifRStr = "M_{4,4}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifRStar1[i].Val2.Key.GetStr() + "," + MotifRStar1[i].Val2.Dat.GetStr() + "," + MotifRStar1[i].Val3.Dat.GetStr() + ")+(" + MotifRStar2[j].Val2.Key.GetStr() + "," + MotifRStar2[j].Val2.Dat.GetStr() + "," + MotifRStar2[j].Val3.Dat.GetStr() + ")+(" + MotifRStar3[l].Val2.Key.GetStr() + "," + MotifRStar3[l].Val2.Dat.GetStr() + "," + MotifRStar3[l].Val3.Dat.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                //pre
                else if (MotifRStar1[i].Val3.Key == 0 && MotifRStar2[j].Val3.Key == 0 && MotifRStar3[l].Val3.Key == 0 &&
                    MotifRStar3[l].Val1.Val2 == 0 && MotifRStar3[l].Val1.Val3 == 1 && MotifRStar3[l].Val1.Val4 == 0 &&
                    MotifRStar2[j].Val1.Val2 == 0 && MotifRStar2[j].Val1.Val3 == 1 &&
                    MotifRStar1[i].Val1.Val2 == 0 &&
                    MotifRStar1[i].Val3.Dat < MotifRStar2[j].Val3.Dat && MotifRStar2[j].Val3.Dat < MotifRStar3[l].Val3.Dat)
                {
                    int index = (5 - 1) * 6 + 3 - 1;
                    MotifRStr = "M_{5,3}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifRStar1[i].Val2.Key.GetStr() + "," + MotifRStar1[i].Val2.Dat.GetStr() + "," + MotifRStar1[i].Val3.Dat.GetStr() + ")+(" + MotifRStar2[j].Val2.Key.GetStr() + "," + MotifRStar2[j].Val2.Dat.GetStr() + "," + MotifRStar2[j].Val3.Dat.GetStr() + ")+(" + MotifRStar3[l].Val2.Key.GetStr() + "," + MotifRStar3[l].Val2.Dat.GetStr() + "," + MotifRStar3[l].Val3.Dat.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifRStar1[i].Val3.Key == 0 && MotifRStar2[j].Val3.Key == 0 && MotifRStar3[l].Val3.Key == 0 &&
                    MotifRStar3[l].Val1.Val2 == 0 && MotifRStar3[l].Val1.Val3 == 1 && MotifRStar3[l].Val1.Val4 == 1 &&
                    MotifRStar2[j].Val1.Val2 == 0 && MotifRStar2[j].Val1.Val3 == 1 &&
                    MotifRStar1[i].Val1.Val2 == 0 &&
                    MotifRStar1[i].Val3.Dat < MotifRStar2[j].Val3.Dat && MotifRStar2[j].Val3.Dat < MotifRStar3[l].Val3.Dat)
                {
                    int index = (5 - 1) * 6 + 4 - 1;
                    MotifRStr = "M_{5,4}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifRStar1[i].Val2.Key.GetStr() + "," + MotifRStar1[i].Val2.Dat.GetStr() + "," + MotifRStar1[i].Val3.Dat.GetStr() + ")+(" + MotifRStar2[j].Val2.Key.GetStr() + "," + MotifRStar2[j].Val2.Dat.GetStr() + "," + MotifRStar2[j].Val3.Dat.GetStr() + ")+(" + MotifRStar3[l].Val2.Key.GetStr() + "," + MotifRStar3[l].Val2.Dat.GetStr() + "," + MotifRStar3[l].Val3.Dat.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifRStar1[i].Val3.Key == 0 && MotifRStar2[j].Val3.Key == 0 && MotifRStar3[l].Val3.Key == 0 &&
                    MotifRStar3[l].Val1.Val2 == 1 && MotifRStar3[l].Val1.Val3 == 0 && MotifRStar3[l].Val1.Val4 == 0 &&
                    MotifRStar2[j].Val1.Val2 == 1 && MotifRStar2[j].Val1.Val3 == 0 &&
                    MotifRStar1[i].Val1.Val2 == 1 &&
                    MotifRStar1[i].Val3.Dat < MotifRStar2[j].Val3.Dat && MotifRStar2[j].Val3.Dat < MotifRStar3[l].Val3.Dat)
                {
                    int index = (5 - 1) * 6 + 5 - 1;
                    MotifRStr = "M_{5,5}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifRStar1[i].Val2.Key.GetStr() + "," + MotifRStar1[i].Val2.Dat.GetStr() + "," + MotifRStar1[i].Val3.Dat.GetStr() + ")+(" + MotifRStar2[j].Val2.Key.GetStr() + "," + MotifRStar2[j].Val2.Dat.GetStr() + "," + MotifRStar2[j].Val3.Dat.GetStr() + ")+(" + MotifRStar3[l].Val2.Key.GetStr() + "," + MotifRStar3[l].Val2.Dat.GetStr() + "," + MotifRStar3[l].Val3.Dat.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifRStar1[i].Val3.Key == 0 && MotifRStar2[j].Val3.Key == 0 && MotifRStar3[l].Val3.Key == 0 &&
                    MotifRStar3[l].Val1.Val2 == 1 && MotifRStar3[l].Val1.Val3 == 0 && MotifRStar3[l].Val1.Val4 == 1 &&
                    MotifRStar2[j].Val1.Val2 == 1 && MotifRStar2[j].Val1.Val3 == 0 &&
                    MotifRStar1[i].Val1.Val2 == 1 &&
                    MotifRStar1[i].Val3.Dat < MotifRStar2[j].Val3.Dat && MotifRStar2[j].Val3.Dat < MotifRStar3[l].Val3.Dat)
                {
                    int index = (5 - 1) * 6 + 6 - 1;
                    MotifRStr = "M_{5,6}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifRStar1[i].Val2.Key.GetStr() + "," + MotifRStar1[i].Val2.Dat.GetStr() + "," + MotifRStar1[i].Val3.Dat.GetStr() + ")+(" + MotifRStar2[j].Val2.Key.GetStr() + "," + MotifRStar2[j].Val2.Dat.GetStr() + "," + MotifRStar2[j].Val3.Dat.GetStr() + ")+(" + MotifRStar3[l].Val2.Key.GetStr() + "," + MotifRStar3[l].Val2.Dat.GetStr() + "," + MotifRStar3[l].Val3.Dat.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifRStar1[i].Val3.Key == 0 && MotifRStar2[j].Val3.Key == 0 && MotifRStar3[l].Val3.Key == 0 &&
                    MotifRStar3[l].Val1.Val2 == 0 && MotifRStar3[l].Val1.Val3 == 0 && MotifRStar3[l].Val1.Val4 == 0 &&
                    MotifRStar2[j].Val1.Val2 == 0 && MotifRStar2[j].Val1.Val3 == 0 &&
                    MotifRStar1[i].Val1.Val2 == 0 &&
                    MotifRStar1[i].Val3.Dat < MotifRStar2[j].Val3.Dat && MotifRStar2[j].Val3.Dat < MotifRStar3[l].Val3.Dat)
                {
                    int index = (6 - 1) * 6 + 3 - 1;
                    MotifRStr = "M_{6,3}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifRStar1[i].Val2.Key.GetStr() + "," + MotifRStar1[i].Val2.Dat.GetStr() + "," + MotifRStar1[i].Val3.Dat.GetStr() + ")+(" + MotifRStar2[j].Val2.Key.GetStr() + "," + MotifRStar2[j].Val2.Dat.GetStr() + "," + MotifRStar2[j].Val3.Dat.GetStr() + ")+(" + MotifRStar3[l].Val2.Key.GetStr() + "," + MotifRStar3[l].Val2.Dat.GetStr() + "," + MotifRStar3[l].Val3.Dat.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifRStar1[i].Val3.Key == 0 && MotifRStar2[j].Val3.Key == 0 && MotifRStar3[l].Val3.Key == 0 &&
                    MotifRStar3[l].Val1.Val2 == 0 && MotifRStar3[l].Val1.Val3 == 0 && MotifRStar3[l].Val1.Val4 == 1 &&
                    MotifRStar2[j].Val1.Val2 == 0 && MotifRStar2[j].Val1.Val3 == 0 &&
                    MotifRStar1[i].Val1.Val2 == 0 &&
                    MotifRStar1[i].Val3.Dat < MotifRStar2[j].Val3.Dat && MotifRStar2[j].Val3.Dat < MotifRStar3[l].Val3.Dat)
                {
                    int index = (6 - 1) * 6 + 4 - 1;
                    MotifRStr = "M_{6,4}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifRStar1[i].Val2.Key.GetStr() + "," + MotifRStar1[i].Val2.Dat.GetStr() + "," + MotifRStar1[i].Val3.Dat.GetStr() + ")+(" + MotifRStar2[j].Val2.Key.GetStr() + "," + MotifRStar2[j].Val2.Dat.GetStr() + "," + MotifRStar2[j].Val3.Dat.GetStr() + ")+(" + MotifRStar3[l].Val2.Key.GetStr() + "," + MotifRStar3[l].Val2.Dat.GetStr() + "," + MotifRStar3[l].Val3.Dat.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifRStar1[i].Val3.Key == 0 && MotifRStar2[j].Val3.Key == 0 && MotifRStar3[l].Val3.Key == 0 &&
                    MotifRStar3[l].Val1.Val2 == 1 && MotifRStar3[l].Val1.Val3 == 1 && MotifRStar3[l].Val1.Val4 == 0 &&
                    MotifRStar2[j].Val1.Val2 == 1 && MotifRStar2[j].Val1.Val3 == 1 &&
                    MotifRStar1[i].Val1.Val2 == 1 &&
                    MotifRStar1[i].Val3.Dat < MotifRStar2[j].Val3.Dat && MotifRStar2[j].Val3.Dat < MotifRStar3[l].Val3.Dat)
                {
                    int index = (6 - 1) * 6 + 3 - 1;
                    MotifRStr = "M_{6,5}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifRStar1[i].Val2.Key.GetStr() + "," + MotifRStar1[i].Val2.Dat.GetStr() + "," + MotifRStar1[i].Val3.Dat.GetStr() + ")+(" + MotifRStar2[j].Val2.Key.GetStr() + "," + MotifRStar2[j].Val2.Dat.GetStr() + "," + MotifRStar2[j].Val3.Dat.GetStr() + ")+(" + MotifRStar3[l].Val2.Key.GetStr() + "," + MotifRStar3[l].Val2.Dat.GetStr() + "," + MotifRStar3[l].Val3.Dat.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
                else if (MotifRStar1[i].Val3.Key == 0 && MotifRStar2[j].Val3.Key == 0 && MotifRStar3[l].Val3.Key == 0 &&
                    MotifRStar3[l].Val1.Val2 == 1 && MotifRStar3[l].Val1.Val3 == 1 && MotifRStar3[l].Val1.Val4 == 1 &&
                    MotifRStar2[j].Val1.Val2 == 1 && MotifRStar2[j].Val1.Val3 == 1 &&
                    MotifRStar1[i].Val1.Val2 == 1 &&
                    MotifRStar1[i].Val3.Dat < MotifRStar2[j].Val3.Dat && MotifRStar2[j].Val3.Dat < MotifRStar3[l].Val3.Dat)
                {
                    int index = (6 - 1) * 6 + 6 - 1;
                    MotifRStr = "M_{6,6}:\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                    MotifRStr = "(" + MotifRStar1[i].Val2.Key.GetStr() + "," + MotifRStar1[i].Val2.Dat.GetStr() + "," + MotifRStar1[i].Val3.Dat.GetStr() + ")+(" + MotifRStar2[j].Val2.Key.GetStr() + "," + MotifRStar2[j].Val2.Dat.GetStr() + "," + MotifRStar2[j].Val3.Dat.GetStr() + ")+(" + MotifRStar3[l].Val2.Key.GetStr() + "," + MotifRStar3[l].Val2.Dat.GetStr() + "," + MotifRStar3[l].Val3.Dat.GetStr() + ")" + "\n";
                    if (!(*motifRecordFile)[index].IsIn(MotifRStr)) { (*motifRecordFile)[index].Add(MotifRStr); printf("\n%s", MotifRStr.CStr()); }
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// Generic three temporal edge, three node star and triad counter.
template <typename EdgeData>
void StarTriad3TEdgeCounter<EdgeData>::Count(const TVec<TThree>& events,
                                             const TIntV& timestamps, double delta, TVec<TVec<TStr>>* motifRecordFile) {
  InitializeCounters();
  if (events.Len() != timestamps.Len()) {
    TExcept::Throw("Number of events must match number of timestamps.");
  }
  int start = 0;
  int end = 0;
  int L = timestamps.Len();
  TVec<TRecordStar> MotifRStar;
  for (int j = 0; j < L; j++) {
    double tj = double(timestamps[j]);
    // Adjust counts in pre-window [tj - delta, tj)
    while (start < L && double(timestamps[start]) < tj - delta) {
      PopPre(events[start], timestamps[start], &MotifRStar);
      start++;
    }
    // Adjust counts in post-window (tj, tj + delta]
    while (end < L && double(timestamps[end]) <= tj + delta) {
      PushPos(events[end], timestamps[end], &MotifRStar);
      end++;
    }
    // Move current event off post-window
    PopPos(events[j], timestamps[j], &MotifRStar);
    ProcessCurrent(events[j], timestamps[j], &MotifRStar, motifRecordFile);
    PushPre(events[j], timestamps[j], &MotifRStar);
  }
}

///////////////////////////////////////////////////////////////////////////////
// Methods for the main sub-routine in the fast star counting algorithm.
void ThreeTEdgeStarCounter::InitializeCounters() {
  pre_sum_ = Counter2D(2, 2);
  pos_sum_ = Counter2D(2, 2);
  mid_sum_ = Counter2D(2, 2);
  pre_counts_ = Counter3D(2, 2, 2);
  pos_counts_ = Counter3D(2, 2, 2);
  mid_counts_ = Counter3D(2, 2, 2);
  pre_nodes_ = Counter2D(2, max_nodes_);
  pos_nodes_ = Counter2D(2, max_nodes_);
}

void ThreeTEdgeStarCounter::PopPre(const TThree& event, const TInt& timestamps, TVec<TRecordStar>* MotifRStar) {
//*MotifRStar=TRecordStar <TIntQu, TIntPair, TInt>:<(level,i,j,event), uv, (lable,timestamps)>
    //lable:0=pre,1=mid,2=pos
    //level:0=nodes,1=sum,2=count
    //j:previous dir
    //event:dir
    int nbr = event.Val1;
    int dir = event.Val2;
    TIntPair uv = event.Val3;
    pre_nodes_(dir, nbr) -= 1;
    (*MotifRStar).DelAll(TRecordStar(TIntQu(0, dir, 100, 100), uv, TIntPair(0,timestamps)));
    for (int i = 0; i < 2; i++) {
        int k = pre_sum_(dir, i);
        pre_sum_(dir, i) -= pre_nodes_(i, nbr);
        if (k > pre_sum_(dir, i)) {
            (*MotifRStar).DelAll(TRecordStar(TIntQu(1, i, dir, 100), uv, TIntPair(0, timestamps)));
        }
    }
}

void ThreeTEdgeStarCounter::PopPos(const TThree& event, const TInt& timestamps, TVec<TRecordStar>* MotifRStar) {
    int nbr = event.Val1;
    int dir = event.Val2;
    TIntPair uv = event.Val3;
    pos_nodes_(dir, nbr) -= 1;
    (*MotifRStar).DelAll(TRecordStar(TIntQu(0, dir, 100, 100), uv, TIntPair(2, timestamps)));
    for (int i = 0; i < 2; i++) { 
        int k = pos_sum_(dir, i);
        pos_sum_(dir, i) -= pos_nodes_(i, nbr);
        if (k > pos_sum_(dir, i)) {
            (*MotifRStar).DelAll(TRecordStar(TIntQu(1, i, dir, 100), uv, TIntPair(2, timestamps)));
        }
    }
}

void ThreeTEdgeStarCounter::PushPre(const TThree& event, const TInt& timestamps, TVec<TRecordStar>* MotifRStar) {
    int nbr = event.Val1;
    int dir = event.Val2;
    TIntPair uv = event.Val3;
    for (int i = 0; i < 2; i++) { 
        int k = pre_sum_(i, dir);
        pre_sum_(i, dir) += pre_nodes_(i, nbr);
        if (k < pre_sum_(i, dir) && !(*MotifRStar).IsIn(TRecordStar(TIntQu(1, i, dir, 100), uv, TIntPair(0, timestamps)))) {
            (*MotifRStar).Add(TRecordStar(TIntQu(1, i, dir, 100), uv, TIntPair(0, timestamps)));
        }
    }
    pre_nodes_(dir, nbr) += 1;
    (*MotifRStar).Add(TRecordStar(TIntQu(0, dir, 100, 100), uv, TIntPair(0, timestamps)));
}

void ThreeTEdgeStarCounter::PushPos(const TThree& event, const TInt& timestamps, TVec<TRecordStar>* MotifRStar) {
    int nbr = event.Val1;
    int dir = event.Val2;
    TIntPair uv = event.Val3;
    for (int i = 0; i < 2; i++) { 
        int k = pos_sum_(i, dir);
        pos_sum_(i, dir) += pos_nodes_(i, nbr);
        if (k < pos_sum_(i, dir) && !(*MotifRStar).IsIn(TRecordStar(TIntQu(1, i, dir, 100), uv, TIntPair(2, timestamps)))) {
            (*MotifRStar).Add(TRecordStar(TIntQu(1, i, dir, 100), uv, TIntPair(2, timestamps)));
        }
    }
    pos_nodes_(dir, nbr) += 1;
    (*MotifRStar).Add(TRecordStar(TIntQu(0, dir, 100, 100), uv, TIntPair(2, timestamps)));
}

void ThreeTEdgeStarCounter::ProcessCurrent(const TThree& event, const TInt& timestamps, TVec<TRecordStar>* MotifRStar, TVec<TVec<TStr>>* motifRecordFile) {
    int nbr = event.Val1;
    int dir = event.Val2;
    TIntPair uv = event.Val3;
    int k = 0;
    // Decrement middle sum
    for (int i = 0; i < 2; i++) {
        k = mid_sum_(i, dir);
        mid_sum_(i, dir) -= pre_nodes_(i, nbr);
        if (k > mid_sum_(i, dir)) {
            (*MotifRStar).DelAll(TRecordStar(TIntQu(1, i, dir, 100), uv, TIntPair(1, timestamps)));
        }
    }
    // Update counts
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            k = pre_counts_(i, j, dir);
            pre_counts_(i, j, dir) += pre_sum_(i, j);
            if (k < pre_counts_(i, j, dir) && !(*MotifRStar).IsIn(TRecordStar(TIntQu(2, i, j, dir), uv, TIntPair(0, timestamps)))) {
                (*MotifRStar).Add(TRecordStar(TIntQu(2, i, j, dir), uv, TIntPair(0, timestamps)));
                PrintStarMotif(MotifRStar, motifRecordFile);
            }
            k = pos_counts_(dir, i, j);
            pos_counts_(dir, i, j) += pos_sum_(i, j);
            if (k < pos_counts_(dir, i, j) && !(*MotifRStar).IsIn(TRecordStar(TIntQu(2, i, j, dir), uv, TIntPair(2, timestamps)))) {
                (*MotifRStar).Add(TRecordStar(TIntQu(2, i, j, dir), uv, TIntPair(2, timestamps)));
                PrintStarMotif(MotifRStar, motifRecordFile);
            }
            k = mid_counts_(i, dir, j);
            mid_counts_(i, dir, j) += mid_sum_(i, j);
            if (k < mid_counts_(i, dir, j) && !(*MotifRStar).IsIn(TRecordStar(TIntQu(2, i, j, dir), uv, TIntPair(1, timestamps)))) {
                (*MotifRStar).Add(TRecordStar(TIntQu(2, i, j, dir), uv, TIntPair(1, timestamps)));
                PrintStarMotif(MotifRStar, motifRecordFile);
            }
        }
    }
    // Increment middle sum
    for (int i = 0; i < 2; i++) { 
        k = mid_sum_(dir, i);
        mid_sum_(dir, i) += pos_nodes_(i, nbr);
        if (k < mid_sum_(dir, i) && !(*MotifRStar).IsIn(TRecordStar(TIntQu(1, i, dir, 100), uv, TIntPair(1, timestamps)))) {
            (*MotifRStar).Add(TRecordStar(TIntQu(1, i, dir, 100), uv, TIntPair(1, timestamps)));
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// Methods for the main sub-routine in the fast triangle counting algorithm.
void ThreeTEdgeTriadCounter::InitializeCounters() {
  pre_nodes_ = Counter3D(2, 2, max_nodes_);
  pos_nodes_ = Counter3D(2, 2, max_nodes_);
  pre_sum_ = Counter3D(2, 2, 2);
  pos_sum_ = Counter3D(2, 2, 2);
  mid_sum_ = Counter3D(2, 2, 2);
  triad_counts_ = Counter3D(2, 2, 2);
}

void ThreeTEdgeTriadCounter::PopPre(const TriadEdgeData& event) {
  int nbr = event.nbr;
  int dir = event.dir;
  int u_or_v = event.u_or_v;
  if (!IsEdgeNode(nbr)) {
    pre_nodes_(dir, u_or_v, nbr) -= 1;
    for (int i = 0; i < 2; i++) {
      pre_sum_(u_or_v, dir, i) -= pre_nodes_(i, 1 - u_or_v, nbr);
    }
  }
}

void ThreeTEdgeTriadCounter::PopPos(const TriadEdgeData& event) {
  int nbr = event.nbr;
  int dir = event.dir;
  int u_or_v = event.u_or_v;  
  if (!IsEdgeNode(nbr)) {  
    pos_nodes_(dir, u_or_v, nbr) -= 1;
    for (int i = 0; i < 2; i++) {
      pos_sum_(u_or_v, dir, i) -= pos_nodes_(i, 1 - u_or_v, nbr);
    }
  }
}

void ThreeTEdgeTriadCounter::PushPre(const TriadEdgeData& event) {
  int nbr = event.nbr;
  int dir = event.dir;
  int u_or_v = event.u_or_v;  
  if (!IsEdgeNode(nbr)) {  
    for (int i = 0; i < 2; i++) {
      pre_sum_(1 - u_or_v, i, dir) += pre_nodes_(i, 1 - u_or_v, nbr);
    }
    pre_nodes_(dir, u_or_v, nbr) += 1;
  }
}

void ThreeTEdgeTriadCounter::PushPos(const TriadEdgeData& event) {
  int nbr = event.nbr;
  int dir = event.dir;
  int u_or_v = event.u_or_v;
  if (!IsEdgeNode(nbr)) {
    for (int i = 0; i < 2; i++) {
      pos_sum_(1 - u_or_v, i, dir) += pos_nodes_(i, 1 - u_or_v, nbr);
    }
    pos_nodes_(dir, u_or_v, nbr) += 1;
  }
}

void ThreeTEdgeTriadCounter::ProcessCurrent(const TriadEdgeData& event) {
  int nbr = event.nbr;
  int dir = event.dir;
  int u_or_v = event.u_or_v;  
  // Adjust middle sums
  if (!IsEdgeNode(nbr)) {
    for (int i = 0; i < 2; i++) {
      mid_sum_(1 - u_or_v, i, dir) -= pre_nodes_(i, 1 - u_or_v, nbr);
      mid_sum_(u_or_v, dir, i) += pos_nodes_(i, 1 - u_or_v, nbr);
    }
  }
  // Update counts
  if (IsEdgeNode(nbr)) {
    // Determine if the event edge is u --> v or v --> u
    int u_to_v = 0;
    if (((nbr == node_u_) && dir == 0) || ((nbr == node_v_) && dir == 1)) {
      u_to_v = 1;
    }
    // i --> j, k --> j, i --> k    
    triad_counts_(0, 0, 0) += mid_sum_(u_to_v,     0, 0)
                           +  pos_sum_(u_to_v,     0, 1)
                           +  pre_sum_(1 - u_to_v, 1, 1);
    // i --> j, k --> i, j --> k
    triad_counts_(1, 0, 0) += mid_sum_(u_to_v,     1, 0)
                           +  pos_sum_(1 - u_to_v, 0, 1)
                           +  pre_sum_(1 - u_to_v, 0, 1);
    // i --> j, j --> k, i --> k
    triad_counts_(0, 1, 0) += mid_sum_(1 - u_to_v, 0, 0)
                           +  pos_sum_(u_to_v,     1, 1)
                           +  pre_sum_(1 - u_to_v, 1, 0);
    // i --> j, i --> k, j --> k
    triad_counts_(1, 1, 0) += mid_sum_(1 - u_to_v, 1, 0)
                           +  pos_sum_(1 - u_to_v, 1, 1)
                           +  pre_sum_(1 - u_to_v, 0, 0);
    // i --> j, k --> j, k --> i
    triad_counts_(0, 0, 1) += mid_sum_(u_to_v,     0, 1)
                           +  pos_sum_(u_to_v,     0, 0)
                           +  pre_sum_(u_to_v,     1, 1);
    // i --> j, k --> i, k --> j
    triad_counts_(1, 0, 1) += mid_sum_(u_to_v,     1, 1)
                           +  pos_sum_(1 - u_to_v, 0, 0)
                           +  pre_sum_(u_to_v,     0, 1);
    // i --> j, j --> k, k --> i
    triad_counts_(0, 1, 1) += mid_sum_(1 - u_to_v, 0, 1)
                           +  pos_sum_(u_to_v,     1, 0)
                           +  pre_sum_(u_to_v,     1, 0);
    // i --> j, i --> k, k --> j
    triad_counts_(1, 1, 1) += mid_sum_(1 - u_to_v, 1, 1)
                           +  pos_sum_(1 - u_to_v, 1, 0)
                           +  pre_sum_(u_to_v,     0, 0);
  }
}
