#ifndef _chi_directed_graph_graph_vertex_h
#define _chi_directed_graph_graph_vertex_h

#include "chi_graph.h"

#include <map>

//###################################################################
/**General implementation of a directed-graph vertex.*/
struct chi_graph::GraphVertex
{
  int id;
  void* context;

  std::set<int> us_edge;
  std::set<int> ds_edge;

  std::map<int,double> us_weights;
  std::map<int,double> ds_weights;

  GraphVertex(int in_id, void* in_context) :
    id(in_id),
    context(in_context)
  {}

  GraphVertex(int in_id) :
    id(in_id),
    context(nullptr)
  {}

  GraphVertex(const GraphVertex& in_v)
  {
    this->id = in_v.id;
    this->context = in_v.context;

    us_edge = in_v.us_edge;
    ds_edge = in_v.ds_edge;
  }

  GraphVertex(GraphVertex&& in_v)
  {
    this->id = in_v.id;
    this->context = in_v.context;

    us_edge = in_v.us_edge;
    ds_edge = in_v.ds_edge;

    in_v.context = nullptr;
  }

  GraphVertex& operator=(const GraphVertex& in_v)
  {
    this->id = in_v.id;
    this->context = in_v.context;

    us_edge = in_v.us_edge;
    ds_edge = in_v.ds_edge;

    return *this;
  }

  GraphVertex& operator=(GraphVertex&& in_v)
  {
    this->id = in_v.id;
    this->context = in_v.context;

    us_edge = in_v.us_edge;
    ds_edge = in_v.ds_edge;

    in_v.context = nullptr;

    return *this;
  }

  bool operator==(const GraphVertex& other)
  {
    return other.id == this->id;
  }
};

#endif