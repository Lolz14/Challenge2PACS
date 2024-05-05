#ifndef HH_SPARSE_MATRIX_HH
#define HH_SPARSE_MATRIX_HH

#include <iostream>
#include <cassert>
#include <complex>
#include <algorithm>
#include <ranges>
#include <iterator>
#include <numeric>
#include <execution>
#include <iomanip>
#include "matrix_traits.hpp"


namespace algebra
{

/*!
*Template class for sparse matrix.
* @param T is the type stored: works with integral, floating and complex types
* @param O is the way of storing elements in the matrix, either RowWise (default) either ColumnWise
*/

template <class T, StorageOrder O = StorageOrder::RowWise>      //default: row-major
class Matrix{

public:

   /*!
   * @brief Constructor for the Matrix class, takes number of rows, columns and possibly format.
   * @param m number of rows
   * @param n number of cols
   * @param compressed true if the matrix is compressed
   */
    Matrix(size_t m=0, size_t n=0, bool compressed = false) 
        :   m_m(m), m_n(n), m_nnz(0), m_compressed(compressed){}


    /*!
    * @brief Getter (rows). 
    * @return private member relative to the number of rows, const
    */
    auto m() const {return m_m;};
    /*!
    * @brief Setter (rows). 
    * @return private member relative to the number of rows, non-const
    */
    auto &m() {return m_m;};


    /*!
    * @brief Getter (cols).
    * @return private member relative to the number of cols, const
    */ 
    auto n() const {return m_n;};
    /*!
    * @brief Setter (cols). 
    * @return private member relative to the number of cols, non-const
    */
    auto &n() {return m_n;};

    
    /*!
    * @brief Getter (non-zero elems).
    * @return number of non zero elements, private
    */
    auto nnz() const {return m_nnz;};

    
    /*!
    * @brief Getter (Uncompressed Matrix).
    * @return a std::map<key_type,StorageOrder::RowWise> or a std::map<key_type,StorageOrder::ColumnWise,IsLessColWise>, hence the uncompressed matrix type
    *         
    */
    auto matUnc() const {return m_mat_uncomp;};

    
    /*!
    * @brief Getter (Compressed Matrix).
    * @return an std::vector<T>, representing the compressed matrix 
    */
    auto matCom() const {return m_val_comp;};

    
    /*!
    * @brief Getter (inner indices).
    * @return a std::vector<std::size_t>, representing the inner indices
    */
    auto innerCom() const {return m_inner;};

    
    /*!
    * @brief Getter (outer indices).
    * @return a std::vector<std::size_t>, representing the outer indices
    */
    auto outerCom() const {return m_outer;};

    
    /*!
    * @brief Getter (state).
    * @return boolean indicating the state of the matrix, const
    */
    bool is_compressed() const {return m_compressed;};
    /*!
    * @brief Setter (state). 
    * @return boolean indicating the state of the matrix, non-const
    */
    bool &is_compressed() {return m_compressed;};


    /*!
    *   @brief Function that clears the buffers used to store the values of the matrix, in both states of the matrix.
    *   @return VOID
    */
    inline void clear_buffer(){
        if (this->is_compressed())
        {
            this->m_mat_uncomp.clear();
        }
        else
        {
            this->m_val_comp.clear();
            this->m_outer.clear();
            this->m_inner.clear();
        }
    }


    /*!
    * @brief Resizes the matrix and assigns new dimensions. 
    * @param m_new: number of rows
    * @param n_new: number of cols
    */
    inline void resize(const size_t m_new, const size_t n_new){
        this->clear_buffer();
        this->m() = m_new;
        this->n() = n_new;
    }


    /*!
    * @brief Checks if a certain row is present in the matrix.
    * @param idx: index that has to be checked
    * @return boolean equal to true if there is at least a non-zero element in the idx-th row; false otherwise. Exception Handling in case of incorrect idx inputed
    */
    bool check_presence_row(const std::size_t &idx) const;


    /*!
    * @brief Checks if a certain column is present in the matrix
    * @param idx: index that has to be checked
    * @return boolean equal to true if there is at least a non-zero element in the idx-th column; false otherwise. Exception Handling in case of incorrect idx inputed
    */
    bool check_presence_col(const std::size_t &idx) const;

    
/*!
 * @brief Checking if element(i,j) is present (non-null).
 *
 * This function checks whether the element at the specified row and column index
 * is non-null in the sparse matrix. It aborts if the provided indices exceed
 * the dimensions of the matrix.
 *
 * @param i The row index of the element to be checked.
 * @param j The column index of the element to be checked.
 * @return true if the element at (i,j) is non-null, false otherwise.
 *         It aborts if the provided indices exceed the dimensions of the matrix.
 */

    bool check_presence(const std::size_t &i, const std::size_t &j) const;

    
/*!
 * @brief Extracting the idx-th row.
 *
 * This function extracts the idx-th row from the sparse matrix. It returns a pair
 * containing two vectors:
 *  - The first vector contains the indices of the columns of the idx-th row that
 *    contain non-zero (nnz) elements.
 *  - The second vector contains the values of the non-zero elements in the idx-th row.
 *
 * @param idx The index of the row to be extracted.
 * @return A pair containing:
 *         - The indices of the columns of the idx-th row that contain non-zero elements.
 *         - The values of the non-zero elements in the idx-th row.
 */

    std::pair<std::vector<std::size_t>,std::vector<T>> get_row(const std::size_t &idx) const;
/*!
 * @brief Extracting the idx-th column.
 *
 * This function extracts the idx-th column from the sparse matrix. It returns a pair
 * containing two vectors:
 *  - The first vector contains the indices of the rows of the idx-th column that
 *    contain non-zero (nnz) elements.
 *  - The second vector contains the values of the non-zero elements in the idx-th column.
 *
 * @param idx The index of the column to be extracted.
 * @return A pair containing:
 *         - The indices of the rows of the idx-th column that contain non-zero elements.
 *         - The values of the non-zero elements in the idx-th column.
 */

    std::pair<std::vector<std::size_t>,std::vector<T>> get_col(const std::size_t &idx) const;

    
 /*!
 * @brief Compresses the matrix.
 *
 * This function compresses the matrix, transitioning from an uncompressed state to a compressed state.
 * It clears the container for the uncompressed state.
 */

    void compress();        


/*!
 * @brief Uncompresses the matrix.
 *
 * This function uncompresses the matrix, transitioning from a compressed state to an uncompressed state.
 * It clears the containers for the compressed state.
 */
    void uncompress();

    
/*!
 * @brief Reads the matrix from a file in Matrix Market Format (MMF).
 *
 * This function reads the matrix from a file stored in the Matrix Market Format (MMF). 
 * It assumes that the comments are located at the beginning of the file and that the file follows MMF standards.
 * The matrix is read in an uncompressed state and returned accordingly.
 *
 * @param file_name The name of the file containing the matrix in MMF.
 * @return True if the matrix is successfully read, false otherwise.
 */

    bool reader_mmf(const std::string &file_name);  


/*!
 * @brief Get element A(i,j) in a read-only manner.
 *
 * This function retrieves the value of element A(i,j) from the matrix in a read-only manner.
 * If the element is non-zero, its value is returned; otherwise, 0 casted to type T is returned.
 *
 * @param i The row index.
 * @param j The column index.
 * @return The value in position (i,j) if non-zero, otherwise 0 casted to type T.
 */
    T operator()(std::size_t i, std::size_t j) const;


/*!
 * @brief Uploads element A(i,j) with a new value.
 *
 * This function uploads the element A(i,j) with a new value. Adding a zero means removing an element.
 * It can update a null element only if the matrix is uncompressed, and it can remove an element only if uncompressed.
 *
 * @param i The row index.
 * @param j The column index.
 * @param newvalue The new value to be inserted.
 * @return If inserted, returns the new value inserted. If not, raises an error since there is an attempt to modify a null element in compressed format.
 */

    T & operator()(std::size_t i, std::size_t j, const T newvalue);


/*!
 * @brief Matrix-vector product.
 *
 * This function computes the matrix-vector product. Both the matrix and the vector are read-only,
 * passed as const lvalue references. It works with complex elements and aborts if the dimensions are wrong.
 * 
 * @param M The algebra::Matrix<T,O> of dimensions m x n.
 * @param v The std::vector<T> of size n.
 * @return A std::vector<T> of size m containing the matrix-vector product.
 */


friend std::vector<T>
    operator * <>(const Matrix<T,O> &M,const std::vector<T> &v);

/*!
 * @brief Matrix-matrix product or matrix-vector product.
 *
 * This function computes the matrix-matrix product or matrix-vector product, being able to handle both
 * depending on the dimensions of the matrices involved. The storage of the result will be the same as
 * that of the first factor. It can determine only on the fly whether it is doing a matrix-matrix product
 * or actually a matrix-vector one; it will compile all the portions of the code. It aborts if the dimensions
 * are incorrect.
 * 
 * @param M The algebra::Matrix<T,O> of dimensions m1 x n1.
 * @param v The algebra::Matrix<T,O1> of dimensions m2 x n2.
 * @return An algebra::Matrix<T,O> of dimensions m1 x n2 as the result of the product between the two matrices.
 */

template<StorageOrder O1>
friend
Matrix<T, O>
operator * (const Matrix<T, O>& M, const Matrix<T, O1>& v) {
    assert(M.n() == v.m());

    bool vector_prod = (v.n() == 1);

    Matrix<T, O> results(M.m(), v.n());

    if (vector_prod) { // Matrix-Vector Product
        if constexpr (O == StorageOrder::RowWise) {
            for (size_t i = 0; i < M.m(); ++i) {
                T temp = static_cast<T>(0);
                auto row_i = M.get_row(i);

                for (size_t j = 0; j < row_i.first.size(); ++j) {
                    temp += row_i.second[j] * v(row_i.first[j], 0);
                }
                results(i, 0, temp);
            }
        } else { // StorageOrder::ColWise
            std::vector<T> temp(M.m(), static_cast<T>(0));

            for (size_t i = 0; i < M.n(); ++i) {
                auto col_i = M.get_col(i);

                for (size_t j = 0; j < col_i.first.size(); ++j) {
                    temp[col_i.first[j]] += col_i.second[j] * v(i, 0);
                }
            }

            for (size_t i = 0; i < M.m(); ++i) {
                results(i, 0, temp[i]);
            }
        }
    } else { // Matrix-Matrix Product
        for (size_t j = 0; j < v.n(); ++j) {
            auto col_j = v.get_col(j);

            for (size_t i = 0; i < M.m(); ++i) {
                auto row_i = M.get_row(i);

                std::vector<size_t> idx;
                idx.reserve(row_i.first.size());

                for (auto el : row_i.first) {
                    if (std::find(col_j.first.begin(), col_j.first.end(), el) != col_j.first.end()) {
                        idx.push_back(el);
                    }
                }

                if (!idx.empty()) {
                    T res = std::transform_reduce(
                        std::execution::par,
                        idx.begin(),
                        idx.end(),
                        static_cast<T>(0),
                        std::plus<T>(),
                        [&row_i, &col_j](auto el) {
                            auto idx_r = std::distance(row_i.first.begin(), std::find(row_i.first.begin(), row_i.first.end(), el));
                            auto idx_c = std::distance(col_j.first.begin(), std::find(col_j.first.begin(), col_j.first.end(), el));
                            return row_i.second[idx_r] * col_j.second[idx_c];
                        }
                    );

                    results(i, j, res);
                }
            }
        }
    }
    return results;
}

    
/*!
 * @brief Stream operator to visualize an algebra::Matrix<T,O>.
 *
 * This function overloads the stream operator to visualize an algebra::Matrix<T,O>.
 * If the matrix is uncompressed, it stores the matrix in the stream object in order to visualize it
 * as a classical table filled with zeros. If compressed, it stores the matrix in the stream object
 * in order to visualize it as three arrays according to CSR/CSC format. It is good for visualizing
 * small matrices.
 * 
 * @param str An lvalue reference to an std::ostream object.
 * @param M The algebra::Matrix<T,O> to be visualized after being stored in the stream object.
 * @return The reference to @param str with the matrix stored.
 */

    friend 
    std::ostream &
    operator << <>(std::ostream & str, Matrix<T,O> const &M);
    

/*!
 * @brief Norm evaluation of the matrix.
 *
 * This function evaluates the norm of the matrix. It is a template function that takes a template parameter N,
 * which indicates an enum for the requested norm. It then performs tag dispatching to the correct function
 * to evaluate the norm as requested.
 * 
 * @tparam N The template parameter that indicates an enum for the requested norm.
 * @return Tag dispatching to the correct function to evaluate the norm as requested.
 */

    template<NormType N>
    auto norm() const{ return norm(NT<N>{});};  
    

private:
    /* Rows.*/
    size_t m_m;
    /* Cols.*/                              
    size_t m_n;
    /* Non-zero element s*/                              
    size_t m_nnz;
    /* Boolean indicating state of the matrix */                            
    bool m_compressed;                       

    /* Uncompressed format container*/
    MatrixUncompressed<T,O> m_mat_uncomp;    

    /* Compressed format container*/
    std::vector<T> m_val_comp;    
    /* Inner indices container */
        
    std::vector<std::size_t> m_inner;        
   /* Outer indices container */
    std::vector<std::size_t> m_outer;        

    /*!
    *Norm One of the matrix.
    * @return Norm One
    */
    auto norm(NT<NormType::One>) const;
    /*!
    *Infinity Norm of the matrix.
    * @return Infinity Norm
    */      
    auto norm(NT<NormType::Infinity>) const;
    /*!
    *Frobenius Norm of the matrix.
    * @return Frobenius Norm
    */
    auto norm(NT<NormType::Frobenius>) const;
};

template <class T, StorageOrder O = StorageOrder::RowWise> 

 std::vector<T> operator*(const Matrix<T, O>& M, const std::vector<T>& v) {
    assert(v.size() == M.n());

    std::vector<T> results(M.m());

    if constexpr (O == StorageOrder::RowWise) {
        for (size_t i = 0; i < M.m(); ++i) {
            auto row_i = M.get_row(i);
            results[i] = std::transform_reduce(std::execution::par,
                                               row_i.first.begin(), row_i.first.end(),
                                               row_i.second.begin(),
                                               static_cast<T>(0),
                                               std::plus<>(),
                                               [&](size_t j, T) {
                                                   if (j < row_i.first.size()) {
                                                       return row_i.second[j] * v[row_i.first[j]];
                                                   } else {
                                                       return static_cast<T>(0);
                                                   }
                                               });
        }
    } else { // COLWISE
        for (size_t i = 0; i < M.n(); ++i) {
            auto col_i = M.get_col(i);
            std::transform(std::execution::par_unseq,
                           col_i.first.begin(), col_i.first.end(),
                           col_i.second.begin(),
                           results.begin(),
                           [&](size_t j, T val) {
                               if (j < results.size()) {
                                   return results[j] + val * v[i];
                               } else {
                                   return val * v[i];
                               }
                           });
        }
    }

    return results;
}
template <class T, StorageOrder O = StorageOrder::RowWise> 
std::ostream& operator<<(std::ostream& str, const Matrix<T, O>& M) {
    if (!M.is_compressed()) {
        for (size_t i = 0; i < M.m(); ++i) {
            for (size_t j = 0; j < M.n(); ++j) {
                str << std::setw(10) << M(i, j) << ' '; // Adjust the width as needed
            }
            str << '\n';
        }
    } else {
        auto format_and_print = [&str](const auto& vec, const std::string& label) {
            str << label << ":\n";
            for (const auto& val : vec) {
                str << std::setw(10) << val << '\n'; // Adjust the width as needed
            }
        };

        format_and_print(M.m_val_comp, "Values");
        format_and_print(M.m_inner, "Inner");
        format_and_print(M.m_outer, "Outer");
    }
    return str;
}
}

#include "utils.hpp"


#endif     //HH_SPARSE_MATRIX_HH