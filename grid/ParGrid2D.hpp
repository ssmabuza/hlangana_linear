#ifndef ParGrid2D_HPP
#define ParGrid2D_HPP

#include <mpi.h>
#include <parmetis.h>

#include <Teuchos_RCP.hpp>
#include <Tpetra_Core.hpp>

#include "StdAfx.hpp"

#include "Grid2D.hpp"


/**
 * The parallel grid built from the serial Grid2D object
 * 
 * \param [in] grid    - a serial mesh to be partitioned
 * \param [in] ep_comm - the Epetra parallel communicator
 * 
 * \param [out] pargrid - a parallel grid object with facilities for creating maps and target maps for assembly
 */
class ParGrid2D {
  
  
  //! Empty Ctor
  ParGrid2D() {}
  
  
public:
  
  
  //! Ctor from a serial mesh with similar interface
  ParGrid2D(const Teuchos::RCP<Grid2D> &grid, const Teuchos::RCP<Teuchos::MpiComm<int>>& tcomm);
  
  
  //! Dtor 
  ~ParGrid2D();
  

  //! Some get methods for computing externally used mesh info
  //@{

  //! computes the center point of an edge
  Point2D getEdgeCenterPoint(const int edge_ind) const;

  
  //! gets the center of mass in an element
  Point2D getElementCenterOfMass(const int elem_ind) const;
  
  
  //! gets the element vertices 
  const std::vector<Point2D> getElementVertices(const int elem_ind) const;
  
  //! get the edge vertices
  const std::vector<Point2D> getEdgeVertices(const int edge_ind) const;
  
  //! get the center of an edge
  Point2D edgeEvalCenterPoint(const int edge_ind) const;
  
  //! computes the outward normal to an element given the edge index
  Point2D getNormalOfEdge(const int elem_ind, const int local_edge_ind) const;
  
  
  //! returns the area of the element
  double getElementArea(const int elem_ind) const;
  

  //! gives the length of the edges  
  double getEdgeLength(const int edge_ind) const;

  
  //! get the number of points
  int numOfPoints() const;
  
  
  //! get the number of elements
  int numOfElements() const;
  
  //@}
  
  
  //! The main data of the mesh
  //@{
  
  //! get the points of the whole mesh
  const std::vector<Point2D>& points() const;
  
  
  //! get the elements of the mesh on the proc including ghost
  const std::vector<Element2D*>& elems() const;
  
  
  //! get the edges of the mesh including ghosts
  const std::vector<Edge2D*>& edges() const;
  
  
  //! local to global index mapping.
  std::vector<int> loc2glob;
  
  
  //! including ghosts 
  std::vector<int> glob2loc; 

  
  //! all points of the grid (serial)
  std::vector<Point2D> m_points;
  
  
  //! all points of the grid on this process (indexes excluding ghosts)
  std::vector<int> m_points_ind_loc;
  

  //! all points of the grid on this process (indexes including ghosts)
  std::vector<int> m_points_ind;
  
  
  //! all the elements of the grid (including ghosts)
  std::vector<Element2D*> m_elems;
  
  
  //! all the edges of the grid on this processor (including ghosts)
  std::vector<Edge2D*> m_edges;
  
  
  //! neighbors sharing an edge with an element (not sharing nodes) 
  std::vector<std::vector<int>> m_neighbors;

  
  //! node to element mapping (exlcuding ghost nodes)
  std::vector<std::map<int,int>> m_stencil;
  
  
  //! ghost elements
  std::vector<int> ghosts;


  //! number of ghost cells per local cell
  std::vector<int> nghosts; 
  
  
  //! number of local neighbours+1 (the cell itsself) per local cell
  std::vector<int> nlocals; 

  
  //! mapping from global elements to local   
  std::vector<int> glob2locs;

  
  //! mapping from global elements to processes
  std::vector<int> glob2proc;

  
  //! buffer?
  std::vector<int> buf;

  
  //! number of local elements (not including ghosts)
  int _n_elems_local;
  
  
  //! number of local elements (including ghosts)
  int _n_elems_local_incl;
  
  
  //! number of elements in the full serial mesh
  int _n_elems;
  
  
  //! own displacement (global PETSC indices start here)
  int displown; 

  
  //@}
  
private:
  
  //! Some set methods for forming the mesh data structures 
  //@{
  
  //! computes the area of an element
  double computeElementArea(const int elem_ind);
  
  
  //! computes the edge length
  double computeEdgeLength(const int edge_ind);
  
  
  //! sets the length for all edges 
  void setLengthEachEdge();
  
  
  //! computes the area for all elems
  void setAreaEachElement();
  
  
  //! Adds an edge
  int addEdge(const int node1, const int node2);
  
  
  //! sets the edges of elems in process
  void setEdges();

  
  //! gets index of edge (-1 if not present)
  int getIndexOfTheEdge(const int node1, const int node2);
  
  //! parallel communicator
  const Teuchos::RCP<Teuchos::MpiComm<int>> m_comm;

  //@}
  
  
protected:
  
  std::vector<std::vector<SearchEdgeStruct>> m_search_edge;
    
};

#endif /** ParGrid2D_HPP */