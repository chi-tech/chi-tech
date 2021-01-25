#ifndef CHI_MESHCONTINUUM_LOCALCELLHANDLER_H_
#define CHI_MESHCONTINUUM_LOCALCELLHANDLER_H_

#include "ChiMesh/Cell/cell.h"

namespace chi_mesh
{

//##################################################
/**Stores references to global cells to enable an iterator.*/
class LocalCellHandler
{
  friend class MeshContinuum;
public:
  std::vector<chi_mesh::Cell*>& native_cells;
  std::vector<chi_mesh::Cell*>& foreign_cells;

private:
  /**Constructor.*/
  LocalCellHandler(std::vector<chi_mesh::Cell*>& in_native_cells,
                   std::vector<chi_mesh::Cell*>& in_foreign_cells) :
    native_cells(in_native_cells),
    foreign_cells(in_foreign_cells)
  {}

public:
  ~LocalCellHandler()
  {
    for (auto& cell : native_cells) delete cell;
    for (auto& cell : foreign_cells) delete cell;
  }

  chi_mesh::Cell& operator[](uint64_t cell_local_index);


  //##################################### iterator Class Definition
  /**Internal iterator class.*/
  class iterator
  {
  public:
    LocalCellHandler& ref_block;
    size_t      ref_element;

    iterator(LocalCellHandler& in_block, size_t i) :
      ref_block(in_block),
      ref_element(i) {}

    iterator operator++() { iterator i = *this; ref_element++; return i; }

    chi_mesh::Cell& operator*() { return *(ref_block.native_cells[ref_element]); }
    chi_mesh::Cell* operator->() { return ref_block.native_cells[ref_element]; }
    bool operator==(const iterator& rhs) const { return ref_element == rhs.ref_element; }
    bool operator!=(const iterator& rhs) const { return ref_element != rhs.ref_element; }
  };

  iterator begin() {return {*this,0};}

  iterator end(){return {*this, native_cells.size()};}

  size_t size() const {return native_cells.size();}
};

}//namespace chi_mesh

#endif //CHI_MESHCONTINUUM_LOCALCELLHANDLER_H_