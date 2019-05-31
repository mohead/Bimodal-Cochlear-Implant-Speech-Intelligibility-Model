#ifndef GTFB_STRICT_FILTER_HH
#define GTFB_STRICT_FILTER_HH

/**
 * @file includes classes for the gammatone filter as described by Hohmann.
 */

#include <cassert>
#include "complex.hh"

namespace Gtfb {
    /*
     * Namespace contains filter classes that implement Hohmann's Gammatone
     * Filters exactly as described in [Hohmann 2002]
     */
    namespace Strict {
        /**
         * A Class implementing an IIR filter of first order with only one
         * pole and no zero.
         */
        template <class F>
        class Filter1 {
        public:
            /**
             *  The current state of this filter
             */
            Complex<F> state;

            /**
             * The constructor takes a complex coefficient as argument.
             */
            inline
            Filter1()
                : state(0)
            {}

            /**
             * This method filters 1 complex sample.
             * @param coefficient    The complex coefficient of the filter.
             * @param input_sample   The complex input sample.
             * @return               The complex output sample. 
             */
            inline
            const Complex<F> & filter_complex(const Complex<F> & coefficient,
                                              const Complex<F> & input_sample)
            {
                this->state.mult_complex(coefficient);
                this->state.add_complex(input_sample);
                return this->state;
            }
            /**
             * This method filters 1 real input sample.  Output is of course
             * still complex.
             * @param coefficient    The complex coefficient of the filter.
             * @param input_sample   The real input sample.
             * @return               The complex output sample. 
             */
            inline
            const Complex<F> & filter_real(const Complex<F> & coefficient,
                                           F input_sample)
            {
                this->state.mult_complex(coefficient);
                this->state.add_real(input_sample);
                return this->state;
            }
        };

        /**
         * A Class implementing an IIR filter of Nth order with only one pole
         * and no zero.
         */
        template <class F>
        class FilterN{
        public:
            /**
             * The complex coefficient used by all stages of this filter.
             */
            Complex<F> coefficient;

            /**
             * The order of this filter
             */
            const size_t order;

            /**
             * An array of 1st order filters used by this Nth order filter.
             */
            Filter1<F> * const stages;
    
            /**
             * Constructs an Nth order all-pole filter with all poles equal to
             * _coefficient.
             * @param _order         The desired order of the filter.
             * @param _coefficient   The coefficient used in all filter stages.
             */
            inline
            FilterN(size_t _order, const Complex<F> & _coefficient)
                : coefficient(_coefficient),
                  order(_order),
                  stages(new Filter1<F>[order])
            {}
            /**
             * Destructor deletes all filter stages.
             */
            ~FilterN()
            {
                delete [] stages;
            }
            /**
             * Returns the 1st order filter at the given stage after asserting
             * that stage_no is in the valid range.
             * @param stage_no   The stage number of the desired fílter
             *                   (0 <= stage_no < order).
             * @return           The desired filter.
             */
            inline
            Filter1<F> & stage(size_t stage_no)
            {
                assert(stage_no < this->order);
                return stages[stage_no];
            }
            /**
             * This method filters 1 complex sample.
             * @param input_sample   The complex input sample.
             * @return               The complex output sample. 
             */
            inline
            const Complex<F> & filter_complex(const Complex<F> & input_sample)
            {
                const Complex<F> * result = &input_sample;
                for (size_t stage_no = 0; stage_no < this->order; ++stage_no) {
                    Filter1<F> & filter1 = this->stage(stage_no);
                    result = &filter1.filter_complex(this->coefficient, *result);
                }
                return *result;
            }
            /**
             * This method filters 1 real sample.  The output is still complex.
             * @param input_sample   The real input sample.
             * @return               The complex output sample. 
             */
            inline
            const Complex<F> & filter_real(F input_sample)
            {
                const Complex<F> * result =
                    &this->stage(0).filter_real(this->coefficient, input_sample);
                for (size_t stage_no = 1; stage_no < this->order; ++stage_no) {
                    Filter1<F> & filter1 = this->stage(stage_no);
                    result = &filter1.filter_complex(this->coefficient, *result);
                }
                return *result;
            }
        };
    }
}


#ifdef SWIG // Wrapper generator used for unit tests
%{
#include "strict_filter.hh"
%}
%template (Strict_Filter1D) Gtfb::Strict::Filter1<double>;
%template (Strict_Filter1F) Gtfb::Strict::Filter1<float>;
%template (Strict_FilterND) Gtfb::Strict::FilterN<double>;
%template (Strict_FilterNF) Gtfb::Strict::FilterN<float>;
#endif        


// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// compile-command: "make -f Testfile"
// End:

#endif
