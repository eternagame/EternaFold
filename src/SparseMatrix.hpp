//////////////////////////////////////////////////////////////////////
// SparseMatrix.hpp
//////////////////////////////////////////////////////////////////////

#ifndef SPARSEMATRIX_HPP
#define SPARSEMATRIX_HPP

#include "Utilities.hpp"

//////////////////////////////////////////////////////////////////////
// struct SparseMatrixEntry
//////////////////////////////////////////////////////////////////////

template<class T>
struct SparseMatrixEntry
{
    int column;
    T value;

    // constructors and assignment operator
    SparseMatrixEntry();

    // constructor
    SparseMatrixEntry(int column, T value);

    // copy constructor
    SparseMatrixEntry(const SparseMatrixEntry &rhs);

    // assignment operator
    SparseMatrixEntry &operator=(const SparseMatrixEntry &rhs);
};

//////////////////////////////////////////////////////////////////////
// class SparseMatrix
//////////////////////////////////////////////////////////////////////

template<class T>
class SparseMatrix
{
    SparseMatrixEntry<T> *data;
    SparseMatrixEntry<T> **row_ptrs;
    int rows, cols, entries;
    T zero;
    
public:
    enum ConversionType { TRANSPOSE };

    // create an empty sparse matrix without allocating memory
    SparseMatrix();

    // create empty sparse matrix, using row_sizes[] to determine the number of entries allocated for each row
    SparseMatrix(const std::vector<int> &row_sizes, int rows, int cols, T zero);

    // create sparse matrix from a map of (row, column) pairs to values
    SparseMatrix(const std::map<std::pair<int,int>, T> &source, int rows, int cols, T zero);

    // create sparse matrix from dense rectangular array, using existing sparse matrix to determine sparsity pattern
    SparseMatrix(const SparseMatrix<T> &mask, const T *source);

    // create sparse matrix from dense rectangular array
    SparseMatrix(const T *source, int rows, int cols, T zero);

    // create sparse matrix from dense upper triangular array
    SparseMatrix(const T *source, int size, T zero);

    // copy constructor
    SparseMatrix(const SparseMatrix &rhs);

    // copy constructor, with transpose
    SparseMatrix(const SparseMatrix &rhs, ConversionType);

    // destructor
    ~SparseMatrix();

    // assignment operator
    SparseMatrix &operator=(const SparseMatrix &rhs);

    // accessors

    // get number of rows
    int GetNumRows() const;

    // get number of columns
    int GetNumCols() const;

    // get number of non-zero entries
    int GetNumEntries() const;

    // get pointer to beginning of row
    const SparseMatrixEntry<T> *GetRowBegin(int row) const;
    
    // get pointer past end of row
    const SparseMatrixEntry<T> *GetRowEnd(int row) const;

    // get pointer to last element of row
    const SparseMatrixEntry<T> *GetRowRBegin(int row) const;

    // get pointer to element before beginning of row
    const SparseMatrixEntry<T> *GetRowREnd(int row) const;

    // get pointer to beginning of row
    SparseMatrixEntry<T> *GetRowBegin(int row);

    // get pointer past end of row
    SparseMatrixEntry<T> *GetRowEnd(int row);

    // get pointer to last element of row
    SparseMatrixEntry<T> *GetRowRBegin(int row);

    // get pointer to element before beginning of row
    SparseMatrixEntry<T> *GetRowREnd(int row);

    // retrieve an arbitrary element
    T operator()(int row, int col) const;

    // compute sum of all non-zero elements
    T GetSum() const;

    // return rows*cols length vector containing all entries
    std::vector<T> GetUnsparse() const;

    // return rows*cols length vector containing 1 for each non-zero entry
    // and zero elsewhere
    std::vector<T> GetUnsparseMask() const;

    // print a sparse version of the matrix to a file
    void PrintSparse(std::ostream &outfile) const;

    // print a sparse version of the matrix, along with sequence letters, to a file
    void PrintSparseBPSEQ(std::ostream &outfile, const std::string &s) const;
};

#include "SparseMatrix.ipp"

#endif
