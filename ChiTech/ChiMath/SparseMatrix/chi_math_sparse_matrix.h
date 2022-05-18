#ifndef _chi_math_sparse_matrix_h
#define _chi_math_sparse_matrix_h

#include "../chi_math.h"

//###################################################################
/**Sparse matrix utility. This is a basic CSR type sparse matrix
 * which allows efficient matrix storage and multiplication. It is
 * not intended for solving linear systems (use PETSc for that instead).
 * It was originally developed for the transfer matrices of transport
 * cross-sections.*/
class chi_math::SparseMatrix
{
private:
  size_t row_size;   ///< Maximum number of rows for this matrix
  size_t col_size;   ///< Maximum number of columns for this matrix

public:
  /**rowI_indices[i] is a vector indices j for the
   * non-zero columns.*/
  std::vector<std::vector<size_t>> rowI_indices;
  /**rowI_values[i] corresponds to column indices and
   * contains the non-zero value.*/
  std::vector<std::vector<double>> rowI_values;

public:
  SparseMatrix(size_t num_rows, size_t num_cols);
  SparseMatrix(const SparseMatrix& in_matrix);

  size_t NumRows() const {return row_size;}
  size_t NumCols() const {return col_size;}

  void   Insert(size_t i, size_t j, double value);
  void   InsertAdd(size_t i, size_t j, double value);
  double ValueIJ(size_t i, size_t j) const;
  void   SetDiagonal(const std::vector<double>& diag);

  void Compress();

  std::string PrintStr() const;

private:
  void CheckInitialized() const;

public:
  virtual ~SparseMatrix() = default;

public:
  class EntryReference
  {
  public:
    const size_t& row_index;
    const size_t& column_index;
    double&       value;

    EntryReference(const size_t& row_id,
                   const size_t& column_id,
                   double& in_value) :
             row_index(row_id),
             column_index(column_id),
             value(in_value) {}
  };

  class ConstEntryReference
  {
  public:
    const size_t& row_index;
    const size_t& column_index;
    const double& value;

    ConstEntryReference(const size_t& row_id,
                        const size_t& column_id,
                        const double& in_value) :
      row_index(row_id),
      column_index(column_id),
      value(in_value) {}
  };

  class RowIteratorContext
  {
  private:
    const std::vector<size_t>& m_ref_col_ids;
          std::vector<double>& m_ref_col_vals;
    const size_t               m_ref_row;
  public:
    RowIteratorContext(SparseMatrix& matrix, size_t ref_row) :
                      m_ref_col_ids(matrix.rowI_indices[ref_row]),
                      m_ref_col_vals(matrix.rowI_values[ref_row]),
                      m_ref_row(ref_row){}

    class RowIterator
    {
    private:
      typedef RowIterator It;
    private:
      RowIteratorContext& m_context;
      size_t              m_ref_entry;
    public:
      RowIterator(RowIteratorContext& context, size_t ref_entry) :
                  m_context{context}, m_ref_entry{ref_entry} {}

      It operator++()    {It i = *this; m_ref_entry++; return i;}
      It operator++(int) {m_ref_entry++; return *this;}

      EntryReference
      operator*() {return {m_context.m_ref_row,
                           m_context.m_ref_col_ids[m_ref_entry],
                           m_context.m_ref_col_vals[m_ref_entry]};}

      bool operator==(const It& rhs) const {return m_ref_entry == rhs.m_ref_entry;}
      bool operator!=(const It& rhs) const {return m_ref_entry != rhs.m_ref_entry;}
    };

    RowIterator begin() {return {*this, 0};}
    RowIterator end()   {return {*this, m_ref_col_vals.size()};}
  };

  RowIteratorContext Row(size_t row_id); //See .cc file

  class ConstRowIteratorContext
  {
  private:
    const std::vector<size_t>& m_ref_col_ids;
    const std::vector<double>& m_ref_col_vals;
    const size_t               m_ref_row;
  public:
    ConstRowIteratorContext(const SparseMatrix& matrix, size_t ref_row) :
      m_ref_col_ids(matrix.rowI_indices[ref_row]),
      m_ref_col_vals(matrix.rowI_values[ref_row]),
      m_ref_row(ref_row){}

    class ConstRowIterator
    {
    private:
      typedef ConstRowIterator It;
    private:
      const ConstRowIteratorContext& m_context;
      size_t              m_ref_entry;
    public:
      ConstRowIterator(const ConstRowIteratorContext& context, size_t ref_entry) :
        m_context(context), m_ref_entry{ref_entry} {}

      It operator++()    {It i = *this; m_ref_entry++; return i;}
      It operator++(int) {m_ref_entry++; return *this;}

      ConstEntryReference
      operator*() {return {m_context.m_ref_row,
                           m_context.m_ref_col_ids[m_ref_entry],
                           m_context.m_ref_col_vals[m_ref_entry]};}

      bool operator==(const It& rhs) const {return m_ref_entry == rhs.m_ref_entry;}
      bool operator!=(const It& rhs) const {return m_ref_entry != rhs.m_ref_entry;}
    };

    ConstRowIterator begin() const {return {*this, 0};}
    ConstRowIterator end() const   {return {*this, m_ref_col_vals.size()};}
  };

  ConstRowIteratorContext Row(size_t row_id) const;

  /**Iterator to loop over all matrix entries.*/
  class EntriesIterator
  {
  private:
    typedef EntriesIterator EIt;
  private:
    SparseMatrix& m_matrix;
    size_t m_ref_row;
    size_t m_ref_col;
  public:

    explicit EntriesIterator(SparseMatrix& context, size_t row) :
      m_matrix{context},m_ref_row{row},m_ref_col(0)
      {}

    void Advance()
    {
      m_ref_col++;
      if (m_ref_col >= m_matrix.rowI_indices[m_ref_row].size())
      {
        m_ref_row++;
        m_ref_col = 0;
        while ((m_ref_row < m_matrix.row_size) and
               (m_matrix.rowI_indices[m_ref_row].empty()))
          m_ref_row++;
      }
    }

    EIt operator++() {EIt i = *this; Advance(); return i;}
    EIt operator++(int) {Advance(); return *this;}

    EntryReference operator*()
    {
      return {m_ref_row,
              m_matrix.rowI_indices[m_ref_row][m_ref_col],
              m_matrix.rowI_values[m_ref_row][m_ref_col]};
    }
    bool operator==(const EIt& rhs) const
    { return (m_ref_row == rhs.m_ref_row) and
             (m_ref_col == rhs.m_ref_col); }
    bool operator!=(const EIt& rhs) const
    { return (m_ref_row != rhs.m_ref_row) or
             (m_ref_col != rhs.m_ref_col); }
  };

  EntriesIterator begin();
  EntriesIterator end();
};


#endif