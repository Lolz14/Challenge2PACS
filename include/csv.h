#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

class csv;
/*!
* Forward definition of csv class
* Inline declaration for Row_escape and flush methods
*/
inline static csv& row_escape(csv& file);
inline static csv& flush(csv& file);

class csv
{
    /*!
    *Filename
    */
    std::ofstream m_fs;
    /*!
    *Boolean to check if the row that is being written is the header
    */
    bool header_bool;
    /*!
    *Separator (default=;)
    */
    const std::string m_separator;
    /*!
    *Escape character
    */
    const std::string m_escape;
    /*!
    *Special character
    */
    const std::string m_special;
public:
    /*! Constructor*/
    csv(const std::string filename, const std::string separator = ";")
        : m_fs()
        , header_bool(true)
        , m_separator(separator)
        , m_escape("\"")
        , m_special("\"")
    {
        /*! Exception Handling for bad initialization*/
        m_fs.exceptions(std::ios::failbit | std::ios::badbit);
        m_fs.open(filename);
    }

    /*!
    * Destructor
    * @brief Writes information contained in the buffer using flush, then closes the file
    */

    ~csv()
    {
    
        flush();
        m_fs.close();
    }
    /*!
    * Flushes, hence writes information contained in the buffer
    */
    void flush()
    {
        m_fs.flush();
    }
    /*!
    * Changes row
    */
    void row_escape()
    {
        m_fs << std::endl;
        header_bool = true;
    }

   /*!
 * @brief Overloaded insertion operator for adding values to a CSV object.
 * 
 * This operator allows adding various types of values to a CSV object, such as strings, characters, and custom types. 
 * 
 * @param val A function pointer or function object that accepts a CSV object and returns a modified CSV object.
 * @return A reference to the modified CSV object after applying the specified operation.
 */
csv& operator << ( csv& (* val)(csv&))
{
    return val(*this);
}

/*!
 * @brief Overloaded insertion operator for adding C-style strings to a CSV object.
 * 
 * This operator allows adding C-style strings (const char*) to a CSV object. 
 * The string is first escaped to handle special characters before being added to the CSV.
 * 
 * @param val The C-style string to be added to the CSV object.
 * @return A reference to the modified CSV object after adding the string.
 */
csv& operator << (const char * val)
{
    return write(escape(val));
}

/*!
 * @brief Overloaded insertion operator for adding std::string objects to a CSV object.
 * 
 * This operator allows adding std::string objects to a CSV object. 
 * The string is first escaped to handle special characters before being added to the CSV.
 * 
 * @param val The std::string object to be added to the CSV object.
 * @return A reference to the modified CSV object after adding the string.
 */
csv& operator << (const std::string & val)
{
    return write(escape(val));
}

/*!
 * @brief Template overloaded insertion operator for adding custom types to a CSV object.
 * 
 * This operator template allows adding custom types to a CSV object. 
 * 
 * @tparam T The type of value to be added to the CSV object.
 * @param val The value of type T to be added to the CSV object.
 * @return A reference to the modified CSV object after adding the value.
 */
template<typename T>
csv& operator << (const T& val)
{
    return write(val);
}

private:



/*!
 * @brief Writes a value to the CSV object, separating it with a separator if not the first element in the row.
 * 
 * This function writes a value of type T to the CSV object. If it's not the first element in the row,
 * it appends the value to the row with a separator. It also updates the header boolean to indicate whether
 * the header row has been written or not.
 * 
 * @tparam T The type of value to be written to the CSV object.
 * @param val The value of type T to be written to the CSV object.
 * @return A reference to the modified CSV object after writing the value.
 */
    template<typename T>
    csv& write (const T& val)
    {
        if (!header_bool)
        {
            m_fs << m_separator;
        }
        else
        {
            header_bool = false;
        }
        m_fs << val;
        return *this;
    }
std::string escape(const std::string & val){
        std::ostringstream result;
        result << '"';
        std::string::size_type to, from = 0u, len = val.length();
        while (from < len &&
                std::string::npos != (to = val.find_first_of(m_special, from)))
        {
            result << val.substr(from, to - from) << m_escape << val[to];
            from = to + 1;
        }
        result << val.substr(from) << '"';
        return result.str();
    }};


inline static csv& row_escape(csv& file)
{
    file.row_escape();
    return file;
}

inline static csv& flush(csv& file)
{
    file.flush();
    return file;
}
