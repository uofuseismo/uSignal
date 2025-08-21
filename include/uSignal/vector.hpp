#ifndef USIGNAL_VECTOR_HPP
#define USIGNAL_VECTOR_HPP
#include <vector>
#include <memory>
#include <boost/align.hpp>
namespace USignal
{
template<class T = double>
/// @brief Defines a vector for use in USignal.  This a bit more general
///        than standard C++ vectors and provides simple functions.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class Vector
{
private:
    using DataTypeT = std::vector<T, boost::alignment::aligned_allocator<T, 64>>;
public:
    using iterator = typename DataTypeT::iterator;
    using const_iterator = typename DataTypeT::const_iterator;
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    Vector();
    /// @brief Copy constructor.
    Vector(const Vector &v);
    /// @brief Move constructor.
    Vector(Vector &&v) noexcept;
    /// @brief Constructs from a given vector.
    explicit Vector(const std::vector<T> &vector); 
    /// @brief Constructs a vector of a given size.
    explicit Vector(size_t n);
    /// @brief Constructs a vector of a given size and fills with default values.
    Vector(size_t n, const T value);
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment.
    Vector& operator=(const Vector &v);
    /// @brief Move assignment.
    Vector& operator=(Vector &&v) noexcept;
    /// @}

    /// @brief True indicates that the vector is empty.
    [[nodiscard]] bool empty() const noexcept;
    /// @brief Resizes the vector to the given size.
    void resize(size_t n);
    /// @brief Resizes the vector to teh given size and fills with the value.
    void resize(size_t n, T value);
    /// @brief Reserves space for the vector.
    void reserve(size_t n);
    /// @result The size of the vector.
    [[nodiscard]] size_t size() const noexcept;
    /// @result A pointer to the internal memory array.
    [[nodiscard]] T *data() noexcept;
    /// @result A pointer to the internal memory array.
    [[nodiscard]] const T *data() const noexcept;
    /// @brief Removes the last element in the vector.
    void pop_back();
    /// @brief Adds a new element at the end of the vector.
    void push_back(T value);

    static int getAlignment() noexcept;

    iterator begin();
    iterator end();
    const_iterator cbegin() const;
    const_iterator cend() const;

    [[nodiscard]] T& operator[](size_t index);
    [[nodiscard]] T& operator[](size_t index) const;
    [[nodiscard]] T& at(size_t index);
    [[nodiscard]] T& at(size_t index) const;
 
    /// @name Destructors
    /// @{

    /// @brief Releases memory.
    void clear() noexcept;
    /// @brief Destructor. 
    ~Vector();
    /// @}
private:
    class VectorImpl;
    std::unique_ptr<VectorImpl> pImpl;
};
/// @result a + \textbf{x}
template<typename T> Vector<T> operator+(T a, const Vector<T> &x);
/// @result \textbf{x} + a
template<typename T> Vector<T> operator+(const Vector<T> &x, T a);
/// @result a \textbf{x}
template<typename T> Vector<T> operator*(T a, const Vector<T> &x);
/// @result \textbf{x} a
template<typename T> Vector<T> operator*(const Vector<T> &x, T a); 
/// @result \textbf{x}/a
template<typename T> Vector<T> operator/(const Vector<T> &x, T a);
}
#endif
