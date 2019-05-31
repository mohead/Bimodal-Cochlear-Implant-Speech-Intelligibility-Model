#ifndef GTFB_COMPLEX_HH
#define GTFB_COMPLEX_HH

#include <iostream>

namespace Gtfb {
    /**
     * A templated class for complex arithmetics. The template parameter can be
     * either double or float
     */
    template <class F> // FLOAT class can be double or float
    class Complex {
    public:
        /** Real part of the complex value. */
        F real;

        /** Imaginary part of the complex value. */
        F imag;

        /**
         * Constructor creates a new Complex.
         * @param _real The real part of the new complex value.
         * @param _imag The imaginary part of the new complex value.
         */
        inline
        Complex(F _real = 0, F _imag = 0)
            : real(_real), imag(_imag)
        {}

        /**
         * Adds another complex value to this one, thereby modifying the
         * reciever.
         * @param other   The complex value to add.
         * @return        A reference to the now modified receiver.
         */
        inline
        Complex<F> & add_complex(const Complex<F> & other)
        {
            this->real += other.real;
            this->imag += other.imag;
            return *this;
        }
        /**
         * Adds a real value to this complex value; modifies the receiver.
         * @param other_real   The real value to add.
         * @return             A reference to the now modified receiver.
         */
        inline
        Complex<F> & add_real(F other_real)
        {
            this->real += other_real;
            return *this;
        }
        /**
         * Multiplies this complex with another complex value; modifies the
         * receiver.
         * @param other   The complex value to multiply with.
         * @return        A reference to the now modified receiver.
         */
        inline
        Complex<F> & mult_complex(const Complex<F> & other)
        {
            F tmp_real = this->real * other.real - this->imag * other.imag;
            this->imag = other.real * this->imag + this->real * other.imag;
            this->real = tmp_real;
            return *this;
        }
        /**
         * Multiplies this complex with a real scalar; modifies the receiver.
         * @param other_real   The real value to multiply with.
         * @return             A reference to the now modified receiver.
         */
        inline
        Complex<F> & mult_real(F other_real)
        {
            this->real *= other_real;
            this->imag *= other_real;
            return *this;
        }
        /**
         * Checks for equality.
         * @param other   The complex value to compare against.
         * @return        True if both complex values are equal, else false.
         */
        inline
        bool equal(const Complex<F> & other)
        {
            return (this->real == other.real) && (this->imag == other.imag);
        }
        /**
         * Compute the square of the absolute value of this complex value.
         * @return The square of the absolute value of this complex value.
         */
        inline
        F abs2(void) const
        {
            return this->real * this->real + this->imag * this->imag;
        }
        /**
         * Compute the absolute value of this complex value.
         * @return The absolute value of this complex value.
         */
        inline
        F abs(void) const
        {
            return static_cast<F>(sqrt(static_cast<double>(this->abs2())));
        }
        /**
         * Create a new Complex with the same value as this.
         * @return   A newly created instance of Complex with the same value as
         *           the receiver.
         */
        inline
        Complex<F> clone(void) const
        {
            return Complex<F>(this->real, this->imag);
        }

        /**
         * Divide by complex divisor.  Modify receiver.
         * @param other   Complex divisor.
         * @return        A reference to the mow modified receiver:
         *                The Complex quotient.
         */
        Complex<F> & div_complex(const Complex<F> & other)
        {
            const F tmp_abs2 = other.abs2();
            F tmp_real = this->real * other.real - this->imag * -other.imag;
            this->imag = other.real * this->imag + this->real * -other.imag;
            this->imag /= tmp_abs2;
            this->real = tmp_real / tmp_abs2;
            return *this;
        }
        /**
         * Replace the value of this Complex with its conjugate.
         * @return A reference to this value.
         */
        inline
        Complex<F> & conjugate(void)
        {
            this->imag = -this->imag;
            return *this;
        }
        /**
         * Replace the value of this Complex with its reciprocal.
         * @return A reference to this Complex.
         */
        inline
        Complex<F> & reciprocal(void)
        {
            this->conjugate();
            return this->mult_real(1 / this->abs2());
        }

        /**
         * Print a representation of this Complex number to the given output
         * stream.
         * @param output_stream   The C++ stream to print on.
         * @return                The C++ stream to which this number was
         *                        printed.
         */
        std::ostream & print_on_stream(std::ostream & output_stream) const
        {
            return output_stream << "(" << real << "+" << imag << "i)";
        }
    };

    /**
     * Convenience function to multiply two Complex values and get a new
     * Complex instance as the result. Does not modify any of its arguments.
     * 
     */
    template <class F>
    inline
    Complex<F> operator *(const Complex<F> & factor1,
                          const Complex<F> & factor2)
    {
        return factor1.clone().mult_complex(factor2);
    }
    
    /**
     * Function to make template operator * testable
     */
    template <class F>
    Complex<F> mult_complex_complex(const Complex<F> & factor1,
                                    const Complex<F> & factor2)
    {
        return factor1 * factor2;
    }
}

#ifdef SWIG // Wrapper generator used for unit tests
%{
#include "complex.hh"
%}
%template (ComplexD)        Gtfb::Complex<double>;
%template (ComplexF)        Gtfb::Complex<float>;
%template (mult_complex_complexD)    Gtfb::mult_complex_complex<double>;
%template (mult_complex_complexF)    Gtfb::mult_complex_complex<float>;
#endif

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// compile-command: "make -f Testfile"
// End:

#endif

