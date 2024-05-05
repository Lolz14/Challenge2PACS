#include "matrix.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <exception>
#include <type_traits>

/*! Definition of the functions that were declared in matrix.hpp
* Computing Norms, Reading matrix from file and all the setters and getters.
*/
template<class T, StorageOrder O>
bool algebra::Matrix<T,O>::reader_mmf(const std::string &file_name) {
    std::ifstream file(file_name);

    if (!file)
        throw std::runtime_error("Error opening file: " + file_name);

    // Ignore comments headers, shifting toward the beginning of the matrix definition
    while (file.peek() == '%') file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    size_t num_row, num_col, num_lines;
    file >> num_row >> num_col >> num_lines;

    if (num_row <= 0 || num_col <= 0 || num_lines < 0)
        throw std::runtime_error("Invalid dimensions or number of non-zero elements in the matrix file");

    // Preallocate memory for the matrix to avoid unnecessary reallocations
    this->resize(num_row, num_col);
    this->uncompress(); // Ensure matrix is in uncompressed state

    // Read and store the matrix elements
    for (size_t i = 0; i < num_lines; ++i) {
        size_t row, col;
        T value;
        file >> row >> col >> value;

        if (row < 1 || row > num_row || col < 1 || col > num_col)
            throw std::runtime_error("Invalid indices for an element in the matrix");

        // Store the element in the matrix
        this->operator()(row - 1, col - 1, std::move(value));
    }

    // Check if the actual number of stored elements matches the declared count
    if (this->nnz() != num_lines) {
        throw std::runtime_error("Mismatch in the number of stored elements and the declared count");
    }

    file.close();

    return true;
}

/*! Norms */

template<class T, StorageOrder O>
auto algebra::Matrix<T,O>::norm(NT<NormType::One>) const {
    if (this->m() == 0 || this->n() == 0)
        return static_cast<double>(0);

    std::vector<double> col_sums(m_n);

    // Calculate column sums
    for (size_t col = 0; col < m_n; ++col) {
        auto col_data = this->get_col(col).second;
        col_sums[col] = std::accumulate(col_data.begin(), col_data.end(), 0.0, [](double sum, const T& val) {
            return sum + std::abs(val);
        });
    }

    // Return the maximum column sum
    return *std::max_element(col_sums.begin(), col_sums.end());
}

template<class T, StorageOrder O>
auto algebra::Matrix<T,O>::norm(NT<NormType::Infinity>) const {
    if (this->m() == 0 || this->n() == 0)
        return static_cast<double>(0);

    std::vector<double> row_sums(m_m);

    // Calculate row sums
    for (size_t row = 0; row < m_m; ++row) {
        auto row_data = this->get_row(row).second;
        row_sums[row] = std::accumulate(row_data.begin(), row_data.end(), 0.0, [](double sum, const T& val) {
            return sum + std::abs(val);
        });
    }

    // Return the maximum row sum
    return *std::max_element(row_sums.begin(), row_sums.end());
}

template<class T, StorageOrder O>
auto algebra::Matrix<T,O>::norm(NT<NormType::Frobenius>) const {
    if (this->m() == 0 || this->n() == 0)
        return static_cast<double>(0);

    double sum_squared_abs = 0.0;

    if (this->is_compressed()) {
        // Compressed format
        sum_squared_abs = std::accumulate(m_val_comp.begin(), m_val_comp.end(), 0.0, [](double sum, const T& val) {
            return sum + std::norm(val);
        });
    } else {
        // Uncompressed format
        sum_squared_abs = std::accumulate(m_mat_uncomp.begin(), m_mat_uncomp.end(), 0.0,
                                           [](double sum, const auto& pair) {
                                               return sum + std::norm(pair.second);
                                           });
    }

    return std::sqrt(sum_squared_abs);
}


template<class T, StorageOrder O>
void algebra::Matrix<T,O>::compress() {
    if (!this->is_compressed()) {
        // Clear compressed data
        m_val_comp.clear();
        m_outer.clear();
        m_inner.clear();

        // Reserve space for compressed data
        m_val_comp.reserve(m_nnz);
        m_outer.reserve(m_nnz);

        if constexpr(O == StorageOrder::RowWise) {
            m_inner.resize(m_m + 1);
        } else {
            m_inner.resize(m_n + 1);
        }

        // Convert uncompressed to compressed format
        size_t index = 0;

        for (size_t i = 0; i < m_m; ++i) {
            if constexpr(O == StorageOrder::RowWise) {
                m_inner[i] = index;
            }

            for (const auto& [key, value] : m_mat_uncomp) {
                if (key[0] == i) {
                    m_val_comp.push_back(value);
                    m_outer.push_back(key[1]);
                    ++index;
                }
            }

            if constexpr(O == StorageOrder::ColumnWise) {
                m_inner[i + 1] = index;
            }
        }

        if constexpr(O == StorageOrder::RowWise) {
            m_inner[m_m] = m_val_comp.size();
        }

        this->is_compressed() = true;
        this->clear_buffer();
    }
}



template<class T, StorageOrder O>
void algebra::Matrix<T,O>::uncompress() {
    if (this->is_compressed()) {
        m_mat_uncomp.clear();

        if constexpr(O == StorageOrder::RowWise) {
            for (size_t i = 0; i < m_m; ++i) {
                if (m_inner[i] < m_inner[i + 1]) {
                    auto begin_row_i = m_inner[i];
                    auto end_row_i = m_inner[i + 1];
                    auto range = std::views::iota(begin_row_i, end_row_i);

                    std::for_each(std::execution::par, range.begin(), range.end(), [&](size_t j) {
                        key_type key_j{i, m_outer[j]};
                        m_mat_uncomp.insert({key_j, m_val_comp[j]});
                    });
                }
            }
        } else {
            for (size_t i = 0; i < m_n; ++i) {
                if (m_inner[i] < m_inner[i + 1]) {
                    auto begin_col_i = m_inner[i];
                    auto end_col_i = m_inner[i + 1];
                    auto range = std::views::iota(begin_col_i, end_col_i);

                    std::for_each(std::execution::par, range.begin(), range.end(), [&](size_t j) {
                        key_type key_j{m_outer[j], i};
                        m_mat_uncomp.insert({key_j, m_val_comp[j]});
                    });
                }
            }
        }

        this->is_compressed() = false;
        this->clear_buffer();
    }
}



template<class T, StorageOrder O>
T algebra::Matrix<T,O>::operator()(std::size_t i, std::size_t j) const {   
    if (!this->check_presence(i, j)) {
        return static_cast<T>(0);
    }

    if (!(this->is_compressed())) {
        key_type elem{i, j};
        return m_mat_uncomp.contains(elem) ? m_mat_uncomp.at(elem) : static_cast<T>(0);
    }

    if constexpr(O == StorageOrder::RowWise) {
        auto row_begin = m_outer.cbegin() + m_inner[i];
        auto row_end = m_outer.cbegin() + m_inner[i + 1];
        auto elem = std::find(std::execution::par, row_begin, row_end, j);
        return elem != row_end ? m_val_comp[std::distance(m_outer.cbegin(), elem)] : static_cast<T>(0);
    } else {
        auto col_begin = m_outer.cbegin() + m_inner[j];
        auto col_end = m_outer.cbegin() + m_inner[j + 1];
        auto elem = std::find(std::execution::par, col_begin, col_end, i);
        return elem != col_end ? m_val_comp[std::distance(m_outer.cbegin(), elem)] : static_cast<T>(0);
    }
}


template<class T, StorageOrder O>
T& algebra::Matrix<T, O>::operator()(std::size_t i, std::size_t j, const T newvalue) {
    // Check if the indexes are valid
    if (i >= m_m || j >= m_n) {
        throw std::out_of_range("Index out of range");
    }

    // Handle uncompressed state
    if (!is_compressed()) {
        key_type key_elem{i, j};
        auto it = m_mat_uncomp.find(key_elem);
        if (it != m_mat_uncomp.end()) {
            // Element already exists, modify it
            it->second = newvalue;
        } else {
            // Element doesn't exist, insert it if newvalue is non-zero
            if (newvalue != T{}) {
                m_nnz++;
                m_mat_uncomp.try_emplace(key_elem, newvalue);
            }
        }
        return m_mat_uncomp[key_elem];
    }

    // Handle compressed state
    if (newvalue != T{}) {
        throw std::logic_error("Cannot modify elements or add new non-zero elements in compressed state");
    }

    // Default value for compressed state
    static T zero_value = T{};
    T& result = zero_value;

    if constexpr (O == StorageOrder::RowWise) {
        auto elem = std::find(std::execution::par,
                              m_outer.cbegin() + m_inner[i], m_outer.cbegin() + m_inner[i + 1],
                              j);
        if (elem != m_outer.cbegin() + m_inner[i + 1]) {
            auto index = std::distance(m_outer.cbegin(), elem);
            result = m_val_comp[index];
        }
    } else {
        auto elem = std::find(std::execution::par,
                              m_outer.cbegin() + m_inner[j], m_outer.cbegin() + m_inner[j + 1],
                              i);
        if (elem != m_outer.cbegin() + m_inner[j + 1]) {
            auto index = std::distance(m_outer.cbegin(), elem);
            result = m_val_comp[index];
        }
    }

    return result;
}




template<class T, StorageOrder O>
bool algebra::Matrix<T, O>::check_presence_row(const std::size_t& idx) const {
    assert(idx < m_m); // Check if the passed index is valid

    if constexpr (O == StorageOrder::RowWise) {
        // For RowWise storage, check if the inner index for the row is within bounds
        return (this->is_compressed()) ? (m_inner[idx] < m_inner[idx + 1]) : (m_mat_uncomp.contains({idx, 0}));
    } else {
        // For ColumnWise storage, check if the row index exists in the outer vector
        return (this->is_compressed()) ? (std::find(std::execution::par, m_outer.cbegin(), m_outer.cend(), idx) != m_outer.cend()) : (std::ranges::find_if(m_mat_uncomp, [&idx](const auto& pair) { return pair.first[0] == idx; }) != m_mat_uncomp.cend());
    }
}


/*!
* Checking presence of col idx-th: at least a nnz element in that col.
*/
template<class T, StorageOrder O>
bool algebra::Matrix<T, O>::check_presence_col(const std::size_t& idx) const {
    assert(idx < m_n); // Check if the passed index is valid wrt the dimension of the matrix

    if constexpr (O == StorageOrder::RowWise) {
        // For RowWise storage, check if the column index exists in the outer vector
        return (this->is_compressed()) ? (std::find(std::execution::par, m_outer.cbegin(), m_outer.cend(), idx) != m_outer.cend()) : (std::ranges::find_if(m_mat_uncomp, [&idx](const auto& pair) { return pair.first[1] == idx; }) != m_mat_uncomp.cend());
    } else {
        // For ColumnWise storage, check if the inner index for the column is within bounds
        key_type lb{ 0, idx };
        key_type up{ 0, idx + 1 };
        return (this->is_compressed()) ? (m_inner[idx] < m_inner[idx + 1]) : (!(std::distance(m_mat_uncomp.upper_bound(lb), m_mat_uncomp.lower_bound(up)) == 0) or m_mat_uncomp.contains(lb));
    }
}



template<class T, StorageOrder O>
bool algebra::Matrix<T, O>::check_presence(const std::size_t& i, const std::size_t& j) const {
    assert(i < m_m && j < m_n); // Check if the passed indices are valid wrt the dimension of the matrix

    // Check if there is any element in the specified row and column
    if (!this->check_presence_row(i) || !this->check_presence_col(j)) {
        return false;
    }

    // Check in uncompressed state
    if (!this->is_compressed()) {
        key_type elem { i, j };
        return m_mat_uncomp.contains(elem);
    }

    // Check in compressed state
    if constexpr (O == StorageOrder::RowWise) {
        // Search if the column index exists within the range of row indices
        auto start_search = m_outer.cbegin() + m_inner[i];
        auto sentinel_search = m_outer.cbegin() + m_inner[i + 1];
        return std::binary_search(start_search, sentinel_search, j);
    } else {
        // Search if the row index exists within the range of column indices
        auto start_search = m_outer.begin() + m_inner[j];
        auto sentinel_search = m_outer.begin() + m_inner[j + 1];
        return std::binary_search(start_search, sentinel_search, i);
    }
}


template<class T, StorageOrder O>
std::pair<std::vector<std::size_t>, std::vector<T>> algebra::Matrix<T, O>::get_row(const std::size_t& idx) const {
    std::vector<std::size_t> col_r_idx; // Column indices of the row
    std::vector<T> val_r_idx;           // Values of the row
    std::size_t nr;                     // Number of elements in the row

    if constexpr (O == StorageOrder::RowWise) {
        if (this->is_compressed()) {
            auto start_col = m_outer.cbegin() + m_inner[idx];
            auto sentinel_col = m_outer.cbegin() + m_inner[idx + 1];
            nr = std::distance(start_col, sentinel_col);

            col_r_idx.reserve(nr);
            val_r_idx.reserve(nr);

            std::copy(std::execution::par, start_col, sentinel_col, std::back_inserter(col_r_idx));
            std::copy(std::execution::par, m_val_comp.begin() + m_inner[idx], m_val_comp.begin() + m_inner[idx + 1], std::back_inserter(val_r_idx));
        } else {
            auto lb = key_type{ idx, 0 };
            auto ub = key_type{ idx + 1, 0 };
            auto start_r = m_mat_uncomp.lower_bound(lb);
            auto sentinel_r = m_mat_uncomp.lower_bound(ub);

            nr = std::distance(start_r, sentinel_r);
            col_r_idx.reserve(nr);
            val_r_idx.reserve(nr);

            std::for_each(start_r, sentinel_r, [&col_r_idx, &val_r_idx](auto& p) { col_r_idx.push_back(p.first[1]); val_r_idx.push_back(p.second); });
        }
    } else {
        if (this->is_compressed()) {
            nr = std::count(std::execution::par, m_outer.cbegin(), m_outer.cend(), idx);

            col_r_idx.reserve(nr);
            val_r_idx.reserve(nr);

            std::size_t k = 0;
            std::size_t check_point_row_cols = 0;

            while (val_r_idx.size() < nr) {
                if (m_outer[k] == idx) {
                    val_r_idx.push_back(m_val_comp[k]);

                    auto col_ptr = std::find_if(m_inner.cbegin() + check_point_row_cols, m_inner.cend(), [&k](auto elem) { return elem > k; });
                    col_r_idx.push_back(std::distance(m_inner.cbegin(), col_ptr - 1));
                    check_point_row_cols += col_r_idx.back();
                }
                k++;
            }
        } else {
            MatrixUncompressed<T, O> temp_map;
            std::ranges::copy_if(m_mat_uncomp, std::inserter(temp_map, temp_map.begin()), [&idx](auto& p) { return p.first[0] == idx; });

            nr = temp_map.size();
            col_r_idx.reserve(nr);
            val_r_idx.reserve(nr);

            std::ranges::for_each(temp_map, [&col_r_idx, &val_r_idx](auto& p) { col_r_idx.push_back(p.first[1]); val_r_idx.push_back(p.second); });

            temp_map.clear();
        }
    }
    return std::make_pair(col_r_idx, val_r_idx);
}

template<class T, StorageOrder O>
std::pair<std::vector<std::size_t>, std::vector<T>>
algebra::Matrix<T, O>::get_col(const std::size_t &idx) const {
    std::vector<std::size_t> row_indices;
    std::vector<T> values;

    if constexpr (O == StorageOrder::ColumnWise) {
        if (this->is_compressed()) {
            auto col_start = m_inner[idx];
            auto col_end = m_inner[idx + 1];

            row_indices.reserve(col_end - col_start);
            values.reserve(col_end - col_start);

            for (std::size_t i = col_start; i < col_end; ++i) {
                row_indices.push_back(m_outer[i]);
                values.push_back(m_val_comp[i]);
            }
        } else {
            for (std::size_t i = 0; i < this->m(); ++i) {
                auto it = std::find_if(m_mat_uncomp.begin(), m_mat_uncomp.end(),
                                       [&](const auto& entry) { return entry.first == std::array{idx, i}; });
                if (it != m_mat_uncomp.end()) {
                    row_indices.push_back(it->first[0]);
                    values.push_back(it->second);
                }
            }
        }
    } else { // RowWise
        if (this->is_compressed()) {
            auto row_start = m_inner[idx];
            auto row_end = m_inner[idx + 1];

            row_indices.reserve(row_end - row_start);
            values.reserve(row_end - row_start);

            for (std::size_t i = row_start; i < row_end; ++i) {
                row_indices.push_back(idx);
                values.push_back(m_val_comp[i]);
            }
        } else {
            for (std::size_t i = 0; i < this->n(); ++i) {
                auto it = std::find_if(m_mat_uncomp.begin(), m_mat_uncomp.end(),
                                       [&](const auto& entry) { return entry.first == std::array{i, idx}; });
                if (it != m_mat_uncomp.end()) {
                    row_indices.push_back(it->first[0]);
                    values.push_back(it->second);
                }
            }
        }
    }

    return std::make_pair(row_indices, values);
}
