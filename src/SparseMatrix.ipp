//////////////////////////////////////////////////////////////////////
// SparseMatrix.cpp
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// SparseMatrixEntry::SparseMatrixEntry()
// SparseMatrixEntry::operator=
//
// Constructors and assignment operator.
//////////////////////////////////////////////////////////////////////

template<class T>
inline SparseMatrixEntry<T>::SparseMatrixEntry()
{}

template<class T>
inline SparseMatrixEntry<T>::SparseMatrixEntry(int column, T value) :
    column(column), value(value)
{
    Assert(column >= 0, "Column number must be nonnegative.");
}

template<class T>
inline SparseMatrixEntry<T>::SparseMatrixEntry(const SparseMatrixEntry<T> &rhs) :
    column(rhs.column), value(rhs.value)
{
    Assert(column >= 0, "Column number must be nonnegative.");
}

template<class T>
inline SparseMatrixEntry<T> &SparseMatrixEntry<T>::operator=(const SparseMatrixEntry<T> &rhs)
{
    if (this != &rhs)
    {
        column = rhs.column;
        value = rhs.value;
        Assert(column >= 0, "Column number must be nonnegative.");
    }
    return *this;
}

//////////////////////////////////////////////////////////////////////
// SparseMatrix::SparseMatrix()
//
// Default constructor.
//////////////////////////////////////////////////////////////////////

template<class T>
SparseMatrix<T>::SparseMatrix() :
    data(NULL), row_ptrs(NULL), rows(0), cols(0), entries(0), zero(0)
{}

//////////////////////////////////////////////////////////////////////
// SparseMatrix::SparseMatrix()
//
// Constructor.  Build sparse matrix from a vector of row sizes.
//////////////////////////////////////////////////////////////////////

template<class T>
SparseMatrix<T>::SparseMatrix(const std::vector<int> &row_sizes, int rows, int cols, T zero) :
    rows(rows), cols(cols), entries(0), zero(zero)
{
    Assert(rows > 0, "Number of rows must be positive.");
    Assert(cols > 0, "Number of colums must be positive.");
    Assert(rows == int(row_sizes.size()), "Mismatch in number of rows.");
    
    // count number of entries

    entries = std::accumulate(row_sizes.begin(), row_sizes.end(), 0);
    
    // allocate memory
    
    data = new SparseMatrixEntry<T>[entries];
    row_ptrs = new SparseMatrixEntry<T>*[rows+1];
    
    // build sparse matrix
    
    std::fill(data, data + entries, SparseMatrixEntry<T>(0, zero));
    row_ptrs[0] = data;
    
    for (size_t i = 0; i < row_sizes.size(); i++)
        row_ptrs[i+1] = row_ptrs[i] + row_sizes[i];
}

//////////////////////////////////////////////////////////////////////
// SparseMatrix::SparseMatrix()
//
// Constructor.  Build sparse matrix from a map of (row, column)
// pairs to values.
//////////////////////////////////////////////////////////////////////

template<class T>
SparseMatrix<T>::SparseMatrix(const std::map<std::pair<int,int>, T> &source, int rows, int cols, T zero) :
    rows(rows), cols(cols), entries(0), zero(zero)
{
    Assert(rows > 0, "Number of rows must be positive.");
    Assert(cols > 0, "Number of colums must be positive.");
    
    // count number of entries
    
    entries = int(source.size());
    
    // allocate memory
    
    data = new SparseMatrixEntry<T>[entries];
    row_ptrs = new SparseMatrixEntry<T>*[rows+1];
    
    // build sparse matrix
    
    SparseMatrixEntry<T> *dest_ptr = data;
    
    int current_row = -1;
    for (typename std::map<std::pair<int,int>, T>::const_iterator iter = source.begin(); iter != source.end(); ++iter)
    {
        Assert(0 <= iter->first.first && iter->first.first < rows, "Invalid row.");
        Assert(0 <= iter->first.second && iter->first.second < cols, "Invalid column.");

        // skip until current row
        
        while (current_row < iter->first.first)
        {
            ++current_row;
            row_ptrs[current_row] = dest_ptr;
        }

        // fill entry
        
        dest_ptr->column = iter->first.second;
        dest_ptr->value = iter->second;
        ++dest_ptr;
    }

    // fill in remainder of row pointers
    
    while (current_row < rows)
    {
        ++current_row;
        row_ptrs[current_row] = dest_ptr;
    }
}

//////////////////////////////////////////////////////////////////////
// SparseMatrix::SparseMatrix()
//
// Constructor.  Create sparse matrix from dense rectangular
// array, using existing sparse matrix to determine sparsity
// pattern.
//////////////////////////////////////////////////////////////////////

template<class T>
SparseMatrix<T>::SparseMatrix(const SparseMatrix<T> &mask, const T *source) :
    rows(mask.rows), cols(mask.cols), entries(mask.entries), zero(mask.zero)
{
    // check if source matrix is empty
    
    if (mask.data == NULL)
    {
        data = NULL;
        row_ptrs = NULL;
        return;        
    }

    Assert(rows > 0, "Number of rows must be positive.");
    Assert(cols > 0, "Number of colums must be positive.");
    Assert(entries >= 0, "Number of entries must be nonnegative.");
    
    // allocate memory
        
    data = new SparseMatrixEntry<T>[entries];
    row_ptrs = new SparseMatrixEntry<T>*[rows+1];
    
    // compute proper row offsets
    
    std::copy(mask.data, mask.data + entries, data);
    for (int i = 0; i <= rows; i++)
        row_ptrs[i] = data + (mask.row_ptrs[i] - mask.data);

    // overwrite data using source

    for (int i = 0; i < rows; i++)
    {
        for (SparseMatrixEntry<T> *ptr = row_ptrs[i]; ptr != row_ptrs[i+1]; ++ptr)
            ptr->value = source[i * cols + ptr->column];
    }
}

//////////////////////////////////////////////////////////////////////
// SparseMatrix::SparseMatrix()
//
// Constructor.  Build sparse matrix from a rectangular matrix.
//////////////////////////////////////////////////////////////////////

template<class T>
SparseMatrix<T>::SparseMatrix(const T *source, int rows, int cols, T zero) :
    rows(rows), cols(cols), entries(0), zero(zero)
{
    Assert(rows > 0, "Number of rows must be positive.");
    Assert(cols > 0, "Number of colums must be positive.");
    
    // count number of entries -- since this constructor uses
    // a dense format, it performs an "equality to zero" test in
    // order to decide which elements should be skipped
    
    entries = rows*cols - std::count(source, source + rows*cols, zero);

    // allocate memory
    
    data = new SparseMatrixEntry<T>[entries];
    row_ptrs = new SparseMatrixEntry<T>*[rows+1];

    // build sparse matrix
    
    const T *source_ptr = source;
    SparseMatrixEntry<T> *dest_ptr = data;
    for (int i = 0; i < rows; i++)
    {
        row_ptrs[i] = dest_ptr;
        for (int j = 0; j < cols; j++)
        {

            // store only non-zero entries
            
            if (*source_ptr != zero)
            {
                dest_ptr->column = j;
                dest_ptr->value = *source_ptr;
                ++dest_ptr;
            }
            
            ++source_ptr;
        }
    }
    
    row_ptrs[rows] = data + entries;
}

//////////////////////////////////////////////////////////////////////
// SparseMatrix::SparseMatrix()
//
// Constructor.  Build sparse matrix from a triangular matrix.
//////////////////////////////////////////////////////////////////////

template<class T>
SparseMatrix<T>::SparseMatrix(const T *source, int size, T zero) :
    rows(size), cols(size), entries(0), zero(zero)
{
    
    Assert(size > 0, "Size of triangular matrix must be positive.");
    
    // count number of entries -- since this constructor uses
    // a dense format, it performs an "equality to zero" test in
    // order to decide which elements should be skipped
    
    const T *source_ptr = source;
    for (int i = 0; i < rows; i++)
    {
        for (int j = i; j < cols; j++)
        {
            if (*source_ptr != zero) entries++;
            ++source_ptr;
        }
    }
    
    // allocate memory
    
    data = new SparseMatrixEntry<T>[entries];
    row_ptrs = new SparseMatrixEntry<T>*[rows+1];
    
    // build sparse matrix
    
    source_ptr = source;
    SparseMatrixEntry<T> *dest_ptr = data;
    for (int i = 0; i < rows; i++)
    {
        row_ptrs[i] = dest_ptr;
        for (int j = i; j < cols; j++)
        {

            // store only non-zero entries

            if (*source_ptr != zero)
            {
                dest_ptr->column = j;
                dest_ptr->value = *source_ptr;
                ++dest_ptr;
            }
            ++source_ptr;
        }
    }
    
    row_ptrs[rows] = data + entries;
}

//////////////////////////////////////////////////////////////////////
// SparseMatrix::SparseMatrix()
//
// Constructor.  Build sparse matrix from an existing sparse
// matrix.
//////////////////////////////////////////////////////////////////////

template<class T>
SparseMatrix<T>::SparseMatrix(const SparseMatrix &rhs) :
    rows(rhs.rows), cols(rhs.cols), entries(rhs.entries), zero(rhs.zero)
{

    // check if source matrix is empty
    
    if (rhs.data == NULL)
    {
        data = NULL;
        row_ptrs = NULL;
    }
    else
    {
        Assert(rows > 0, "Number of rows must be positive.");
        Assert(cols > 0, "Number of colums must be positive.");
        Assert(entries >= 0, "Number of entries must be nonnegative.");

        // allocate memory
        
        data = new SparseMatrixEntry<T>[entries];
        row_ptrs = new SparseMatrixEntry<T>*[rows+1];

        // compute proper row offsets
        
        std::copy(rhs.data, rhs.data + entries, data);
        for (int i = 0; i <= rows; i++)
            row_ptrs[i] = data + (rhs.row_ptrs[i] - rhs.data);
    }
}

//////////////////////////////////////////////////////////////////////
// SparseMatrix::SparseMatrix()
//
// Constructor.  Build transpose of sparse matrix from an
// existing sparse matrix.
//////////////////////////////////////////////////////////////////////

template<class T>
SparseMatrix<T>::SparseMatrix(const SparseMatrix &rhs, ConversionType) :
    rows(rhs.cols), cols(rhs.rows), entries(rhs.entries), zero(rhs.zero)
{
    
    Assert(rows > 0, "Number of rows must be positive.");
    Assert(cols > 0, "Number of colums must be positive.");
    Assert(entries >= 0, "Number of entries must be nonnegative.");
    
    data = new SparseMatrixEntry<T>[entries];
    row_ptrs = new SparseMatrixEntry<T>*[rows+1];
    
    // compute number of elements per row
    
    int *row_size = new int[rows+1];
    std::fill(row_size, row_size+rows+1, 0);  
    for (int i = 0; i < rhs.rows; i++)
    {
        for (SparseMatrixEntry<T> *ptr = rhs.row_ptrs[i]; ptr != rhs.row_ptrs[i+1]; ++ptr)
        {
            ++(row_size[ptr->column]);
        }
        
    }
    
    // compute row pointers
    
    row_ptrs[0] = data;
    for (int i = 1; i <= rows; i++)
        row_ptrs[i] = row_ptrs[i-1] + row_size[i-1];
    
    delete [] row_size;
    
    // build sparse matrix
    
    for (int i = 0; i < rhs.rows; i++)
    {
        for (SparseMatrixEntry<T> *ptr = rhs.row_ptrs[i]; ptr != rhs.row_ptrs[i+1]; ++ptr)
        {
            row_ptrs[ptr->column]->column = i;
            row_ptrs[ptr->column]->value = ptr->value;
            ++(row_ptrs[ptr->column]);
        }
    }
    
    // shift row pointers back
    
    for (int i = rows; i >= 1; i--)
        row_ptrs[i] = row_ptrs[i-1];
    row_ptrs[0] = data;  
}

//////////////////////////////////////////////////////////////////////
// SparseMatrix::~SparseMatrix()
//
// Destructor.
//////////////////////////////////////////////////////////////////////

template<class T>
SparseMatrix<T>::~SparseMatrix()
{
    delete [] data;
    delete [] row_ptrs;
}

//////////////////////////////////////////////////////////////////////
// SparseMatrix::SparseMatrix()
//
// Constructor.  Build sparse matrix from an existing sparse
// matrix.
//////////////////////////////////////////////////////////////////////

template<class T>
SparseMatrix<T> &SparseMatrix<T>::operator=(const SparseMatrix &rhs)
{
    if (this != &rhs)
    {
        rows = rhs.rows;
        cols = rhs.cols;
        entries = rhs.entries;
        zero = rhs.zero;
        
        if (rhs.data == NULL)
        {
            data = NULL;
            row_ptrs = NULL;
        }
        else
        {
            Assert(rows > 0, "Number of rows must be positive.");
            Assert(cols > 0, "Number of colums must be positive.");
            Assert(entries >= 0, "Number of entries must be nonnegative.");
            
            data = new SparseMatrixEntry<T>[entries];
            row_ptrs = new SparseMatrixEntry<T>*[rows+1];
            
            std::copy(rhs.data, rhs.data + entries, data);
            for (int i = 0; i <= rows; i++)
                row_ptrs[i] = data + (rhs.row_ptrs[i] - rhs.data);
        }
    }
    
    return *this;
}
 
//////////////////////////////////////////////////////////////////////
// SparseMatrix::GetNumRows()
// SparseMatrix::GetNumCols()
// SparseMatrix::GetNumEntries()
// SparseMatrix::GetRowBegin()
// SparseMatrix::GetRowEnd()
// SparseMatrix::GetRowRBegin()
// SparseMatrix::GetRowREnd()
// SparseMatrix::operator()
// SparseMatrix::GetSum()
//
// Accessors.
//////////////////////////////////////////////////////////////////////

template<class T>
inline int SparseMatrix<T>::GetNumRows() const { return rows; }

template<class T>
inline int SparseMatrix<T>::GetNumCols() const { return cols; }

template<class T>
inline int SparseMatrix<T>::GetNumEntries() const { return entries; }

template<class T>
inline SparseMatrixEntry<T> *SparseMatrix<T>::GetRowBegin(int row)
{
    Assert(row >= 0 && row < rows, "Invalid row.");
    return row_ptrs[row];
}

template<class T>
inline SparseMatrixEntry<T> *SparseMatrix<T>::GetRowEnd(int row)
{
    Assert(row >= 0 && row < rows, "Invalid row.");
    return row_ptrs[row+1];
}

template<class T>
inline SparseMatrixEntry<T> *SparseMatrix<T>::GetRowRBegin(int row)
{
    Assert(row >= 0 && row < rows, "Invalid row.");
    return row_ptrs[row+1]-1;
}

template<class T>
inline SparseMatrixEntry<T> *SparseMatrix<T>::GetRowREnd(int row)
{
    Assert(row >= 0 && row < rows, "Invalid row.");
    return row_ptrs[row]-1;
}

template<class T>
inline const SparseMatrixEntry<T> *SparseMatrix<T>::GetRowBegin(int row) const 
{
    Assert(row >= 0 && row < rows, "Invalid row.");
    return row_ptrs[row];
}

template<class T>
inline const SparseMatrixEntry<T> *SparseMatrix<T>::GetRowEnd(int row) const 
{
    Assert(row >= 0 && row < rows, "Invalid row.");
    return row_ptrs[row+1];
}

template<class T>
inline const SparseMatrixEntry<T> *SparseMatrix<T>::GetRowRBegin(int row) const 
{
    Assert(row >= 0 && row < rows, "Invalid row.");
    return row_ptrs[row+1]-1;
}

template<class T>
inline const SparseMatrixEntry<T> *SparseMatrix<T>::GetRowREnd(int row) const 
{
    Assert(row >= 0 && row < rows, "Invalid row.");
    return row_ptrs[row]-1;
}

template<class T>
inline T SparseMatrix<T>::operator()(int row, int col) const 
{
    Assert(row >= 0 && row < rows, "Invalid row.");
    Assert(col >= 0 && col < cols, "Invalid row.");
    
    // binary search for correct element
    
    SparseMatrixEntry<T> *left = row_ptrs[row];
    SparseMatrixEntry<T> *right = row_ptrs[row+1] - 1;
    while (left <= right)
    {
        SparseMatrixEntry<T> *mid = left + (right - left) / 2;
        if (mid->column > col)
            right = mid-1;
        else if (mid->column < col)
            left = mid+1;
        else
            return mid->value;
    }
    return zero;
}

template<class T>
inline T SparseMatrix<T>::GetSum() const
{
    T sum = zero;
    for (int i = 0; i < entries; i++)
        sum += data[i].value;
    return sum;
}

//////////////////////////////////////////////////////////////////////
// SparseMatrix::GetUnsparse()
//
// Return copy of matrix as a 2D array.
//////////////////////////////////////////////////////////////////////

template<class T>
std::vector<T> SparseMatrix<T>::GetUnsparse() const 
{
    std::vector<T> ret(rows * cols);
    for (int i = 0; i < rows; i++)
    {
        for (const SparseMatrixEntry<T> *p = row_ptrs[i]; p != row_ptrs[i+1]; ++p)
            ret[i * cols + p->column] = p->value;
    }
    return ret;
}

//////////////////////////////////////////////////////////////////////
// SparseMatrix::GetUnsparseMask()
//
// Return 2D array of present positions.
//////////////////////////////////////////////////////////////////////

template<class T>
std::vector<T> SparseMatrix<T>::GetUnsparseMask() const 
{
    std::vector<T> ret(rows * cols);
    for (int i = 0; i < rows; i++)
    {
        for (const SparseMatrixEntry<T> *p = row_ptrs[i]; p != row_ptrs[i+1]; ++p)
            ret[i * cols + p->column] = T(1);
    }
    return ret;
}

//////////////////////////////////////////////////////////////////////
// SparseMatrix::PrintSparse()
//
// Print in sparse matrix format.
//////////////////////////////////////////////////////////////////////

template<class T>
void SparseMatrix<T>::PrintSparse(std::ostream &outfile) const
{
    for (int i = 0; i < rows; i++)
    {
        if (row_ptrs[i] == row_ptrs[i+1]) continue;
        outfile << i;
        for (SparseMatrixEntry<T> *ptr = row_ptrs[i]; ptr != row_ptrs[i+1]; ++ptr)
            outfile << ' ' << ptr->column << ':' << ptr->value;
        outfile << std::endl;
    }
}

//////////////////////////////////////////////////////////////////////
// SparseMatrix::PrintSparseBPSEQ()
//
// Print BPSEQ posteriors in sparse matrix format.
//////////////////////////////////////////////////////////////////////

template<class T>
void SparseMatrix<T>::PrintSparseBPSEQ(std::ostream &outfile, const std::string &s) const 
{
    Assert(int(s.length()) == rows, "Sequence length does not match sparse matrix size.");
    for (int i = 1; i < rows; i++)
    {
        outfile << i << ' ' << s[i];
        for (SparseMatrixEntry<T> *ptr = row_ptrs[i]; ptr != row_ptrs[i+1]; ++ptr)
            outfile << ' ' << ptr->column << ':' << ptr->value;
        outfile << std::endl;
    }
}
