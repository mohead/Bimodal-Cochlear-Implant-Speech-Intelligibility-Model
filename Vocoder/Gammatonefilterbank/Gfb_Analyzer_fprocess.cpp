/* 
 * This file is part of the gammatone filterbank reference implementation
 * described in V. Hohmann's `acta acustica' article.
 *
 * It contains the matlab interface function for the Gfb_Analyzer_fprocess
 * function. See file "README_extension.txt" for a description
 * how to use the C function with matlab.
 * Gfb_Analyzer_fprocess is a Matlab extension in the C programming language
 * with the same functionality as Gfb_Analyzer_process: It performs the
 * filterbank analysis as described in the `acta acustica' article.
 *
 * filename : Gfb_Analyzer_fprocess
 * copyright: Universitaet Oldenburg
 * author   : tp
 * date     : Jan 2002 Sep Oct 2003
 *
 * update   : 
 */

/*-----------------------------------------------------------------------------
 *   Copyright (C) 2002 2003  AG Medizinische Physik,
 *                        Universitaet Oldenburg, Germany
 *                        http://www.physik.uni-oldenburg.de/docs/medi
 *   
 *   Permission to use, copy, and distribute this software/file and its
 *   documentation for any purpose without permission by UNIVERSITAET OLDENBURG
 *   is not granted.
 *   
 *   Permission to use this software for academic purposes is generally
 *   granted.
 *
 *   Permission to modify the software is granted, but not the right to
 *   distribute the modified code.
 *
 *   This software is provided "as is" without expressed or implied warranty.
 *
 *   Author: Tobias Peters (tpeters@uni-oldenburg.de)
 *---------------------------------------------------------------------------*/

/*
 * file: Gfb_Analyzer_fprocess.cpp
 * 
 * Authors: 2001 2002 Tobias Peters <tpeters@uni-oldenburg.de>
 *
 */

#include <string.h>
#include "mex.h"
#include <vector>
#include "strict_filter.hh"

/**
 * This class wraps the input data to the filterbank.
 */
class Filterbank_Input {
public:
  /**
   * The real part of the matlab matrix containing the input data
   */
  const double * const m_real_input;

  /**
   * The imaginary part of the matlab matrix containing the input data
   */
  const double * const m_imag_input;

  /**
   * A flag indicating wether each filterbank channel has the same input
   * or all filterbank channels have different inputs.
   */
  bool m_multiple_inputs;

  /**
   * A flag indicating wether the input to the filterbank is complex
   */
  bool m_complex_input;
  
  /**
   * The number of input channels
   */
  size_t m_input_channels;

  /**
   * The number of samples in each input channel
   */
  size_t m_input_samples;

  Filterbank_Input(const mxArray * input) throw (const char *)
    : m_real_input(mxGetPr(input)),
      m_imag_input(mxGetPi(input)),
      m_multiple_inputs(mxGetM(input) != 1),
      m_complex_input(mxIsComplex(input)),
      m_input_channels(mxGetM(input)),
      m_input_samples(mxGetN(input))
  { 
    if (!mxIsNumeric(input) || mxIsSparse(input) || !mxIsDouble(input)) {
      throw ("Gfb_Analyzer_fprocess: 2nd argument must be a matrix or a "
	     "row vector");
    }
  }
  /**
   * returns the real part of the specified sample number for the desired 
   * channel
   */
  inline
  double real(size_t sample, size_t channel = 0) const throw (const char *) {
    if (m_multiple_inputs == false) {
      channel = 0;
    }
    if (channel >= m_input_channels) {
      throw ("Gfb:Analyzer:fprocess:Filterbank_Input:real: "
	     "channel number too large.");
    }
    if (sample >= m_input_samples) {
      throw ("Gfb:Analyzer:fprocess:Filterbank_Input:real: "
	     "sample number too large.");
    }
    return m_real_input[sample * m_input_channels + channel];
  }
  /** 
   * returns the imaginary part of the specified sample number for the
   * desired channel
   */
  inline
  double imag(size_t sample, size_t channel = 0) const throw (const char *) {
    if (m_multiple_inputs == false) {
      channel = 0;
    }
    if (channel >= m_input_channels) {
      throw ("Gfb:Analyzer:fprocess:Filterbank_Input:imag: "
	     "channel number too large.");
    }
    if (sample >= m_input_samples) {
      throw ("Gfb:Analyzer:fprocess:Filterbank_Input:imag: "
	     "sample number too large.");
    }
    if (m_complex_input == false) {
      return 0.0;
    }
    return m_imag_input[sample * m_input_channels + channel];
  }
};

/**
 * struct Filterbank_Data helps us by bundling the data extracted from
 * the Matlab Gfb_Analyzer object in a single object, with easy to use
 * constructor and destructor functions.
 */
struct Filterbank_Data {
  /**
   * A vector containing the normalization factors of the individual channels.
   */
  double * normalization_factors   ;

  /**
   * two vectors holding the real and the imaginary part, respectively, of
   * each channel's filter coefficients.
   */ 
  double * real_filter_coefficients; double * imag_filter_coefficients;

  /**
   * two vectors holding the real and the imaginary part, respectively, of
   * each channels state. The filter state part corresponding to a specific
   * channel and filter stage is stored at index
   * [channel * filter_order + filter_stage], where filter_order denotes
   * the order of the gammatone filters (usually == 4).
   */   
  double * real_filter_state       ; double * imag_filter_state       ;

  /**
   * The order of the gammatone filters.
   */
  size_t gamma_order;

  /**
   * The number of gammatone filters.
   */
  size_t channels;
};

/**
 * The constructor function for a struct Filterbank_Data. The allocated
 * vectors are large enough to hold the contained arrays for the given number
 * of channels and filter order.
 */
static
struct Filterbank_Data *
Filterbank_Data_new(unsigned channels, unsigned gamma_order)
{
  struct Filterbank_Data * fbd =
    (struct Filterbank_Data *) mxMalloc(1 * sizeof(struct Filterbank_Data));
  fbd->normalization_factors    = (double*) mxMalloc(channels*sizeof(double));
  fbd->real_filter_coefficients = (double*) mxMalloc(channels*sizeof(double));
  fbd->imag_filter_coefficients = (double*) mxMalloc(channels*sizeof(double));
  fbd->real_filter_state        = (double*) mxMalloc(channels*gamma_order*
						     sizeof(double));
  fbd->imag_filter_state        = (double*) mxMalloc(channels*gamma_order*
						     sizeof(double));
  fbd->gamma_order = gamma_order;
  fbd->channels = channels;
  return fbd;
}

/**
 * The destructor of struct Filterbank_Data. It frees all memory
 * occupied by the structure and the contained vectors at once.
 */
static void
Filterbank_Data_delete(struct Filterbank_Data * fbd)
{
  if (fbd == NULL) return;
  mxFree(fbd->normalization_factors);    fbd->normalization_factors    = NULL;
  mxFree(fbd->real_filter_coefficients); fbd->real_filter_coefficients = NULL;
  mxFree(fbd->imag_filter_coefficients); fbd->imag_filter_coefficients = NULL;
  mxFree(fbd->real_filter_state);        fbd->real_filter_state        = NULL;
  mxFree(fbd->imag_filter_state);        fbd->imag_filter_state        = NULL;
  mxFree(fbd);
}

class Filterbank {
public:
  typedef Gtfb::Strict::FilterN<double> FilterN ;
  typedef Gtfb::Complex<double>         Complex ;

  size_t m_channels;
  mxArray ** m_output;
  double * m_real_output;
  double * m_imag_output;
  const mxArray * m_analyzer_in;
  mxArray ** m_analyzer_out;

  Filterbank_Data * m_filterbank_data;
  Filterbank_Input * m_filterbank_input;
  std::vector<FilterN *> m_filter;
  
  ~Filterbank() {
    m_analyzer_out = 0;
    m_output = 0;
    m_real_output = 0;
    m_imag_output = 0;
    m_analyzer_in = 0;
    Filterbank_Data_delete(m_filterbank_data);
    m_filterbank_data = 0;
    delete m_filterbank_input;
    m_filterbank_input = 0;
    for (size_t filter_index = 0;
	 filter_index < m_filter.size();
	 ++filter_index) {
      delete m_filter[filter_index];
      m_filter[filter_index] = 0;
    }
  }

  Filterbank(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    throw (const char *)
    : m_channels(0), m_output(0), m_real_output(0), m_imag_output(0),
      m_analyzer_in(0), m_analyzer_out(0),
      m_filterbank_data(0), m_filterbank_input(0)
  {
    /* Check for proper number of arguments and extract values. */
    if (nrhs != 2) {
      throw "Gfb_Analyzer_fprocess: needs 2 input arguments";
    }
    if (nlhs > 2) {
      throw("Gfb_Analyzer_fprocess: returns one or two parameters");
      return;
    }
    m_output = &plhs[0];
    m_analyzer_out = (nlhs >= 2) ? &plhs[1] : static_cast<mxArray**>(0);

    m_filterbank_input = new Filterbank_Input(prhs[1]);

    m_analyzer_in = prhs[0];
    check_input_analyzer(m_analyzer_in);
    m_channels =
      extract_and_check_channel_no(m_analyzer_in,
				   m_filterbank_input->m_input_channels);
    m_filterbank_data = extract_filterbank_data(m_analyzer_in,
						m_channels);
    /* allocate output matrix */
    *m_output = mxCreateDoubleMatrix(m_channels,
				     m_filterbank_input->m_input_samples,
				     mxCOMPLEX);
    m_real_output = mxGetPr(*m_output);
    m_imag_output = mxGetPi(*m_output);
  }

  void check_input_analyzer(const mxArray * analyzer_in) const
    throw (const char *)
  {
    if (mxGetNumberOfElements(analyzer_in) != 1) {
      throw("Gfb_Analyzer_fprocess: 1st argument must be a single "
	    "Gfb_Analyzer structure");
    }
    /*
     * Check if analyzer_in contains a field `type' with the value
     * "Gfb_Analyzer"
     */
    mxArray * analyzer_in_type = mxGetField(analyzer_in, 0, "type");
    if (analyzer_in_type == 0 || !mxIsChar(analyzer_in_type)) {
      throw("Gfb_Analyzer_fprocess: 1st argument must be a "
	    "Gfb_Analyzer structure");
    }
    char * analyzer_in_typename = mxArrayToString(analyzer_in_type);
    int analyzer_is_Gfb_Analyzer =
      (strcmp(analyzer_in_typename, "Gfb_Analyzer") == 0);
    mxFree(analyzer_in_typename);
    if (!analyzer_is_Gfb_Analyzer) {
      throw("Gfb_Analyzer_fprocess: 1st argument must be a "
	    "Gfb_Analyzer structure");
    }
  }
  size_t extract_and_check_channel_no(const mxArray * analyzer_in,
				      size_t input_channels)
    throw (const char *)
  {
    /* Get the number of gammatone filter channels */
    const mxArray * center_frequencies_hz =
      mxGetField(analyzer_in, 0, "center_frequencies_hz");
    if (center_frequencies_hz == NULL
	|| !mxIsNumeric(center_frequencies_hz)
	|| mxIsComplex(center_frequencies_hz)
	|| mxIsSparse(center_frequencies_hz) 
	|| !mxIsDouble(center_frequencies_hz)
	|| mxGetM(center_frequencies_hz) != 1) {
      throw("Gfb_Analyzer_fprocess: 1st argument must be a Gfb_Analyzer\n"
	    "       structure, and must contain a real column vector\n"
	    "       field named `center_frequencies_hz'.  Use the\n"
	    "       Gfb_Analyzer_new function to construct a valid\n"
	    "       Gfb_Analyzer structure");
    }
    size_t channels = mxGetN(center_frequencies_hz);
    if (input_channels != 1 &&	input_channels != channels) {
      throw("Gfb_Analyzer_fprocess: Number of rows in input data (2nd\n"
	    "       argument) is not 1, and does not match the number\n"
	    "       of bands of the Gfb_Analyzer struct.");
    }
    return channels;
  }

  Filterbank_Data * extract_filterbank_data(const mxArray * analyzer_in,
					    size_t channels)
  {
    const mxArray * analyzer_in_filters =
      mxGetField(analyzer_in, 0, "filters");
    if (analyzer_in_filters == NULL
	|| !mxIsStruct(analyzer_in_filters)
	|| size_t(mxGetNumberOfElements(analyzer_in_filters)) != channels) {
      throw("Gfb_Analyzer_fprocess: 1st argument must be a Gfb_Analyzer\n"
	    "       structure, and must contain a structure field named\n"
	    "       `filters', containing the filter specifications.  Use\n"
	    "       the Gfb_Analyzer_new function to construct a valid\n"
	    "       Gfb_Analyzer structure");
    }

    /* extract the filter order of the gammatone filters */
    mxArray * analyzer_in_gamma_order =
      mxGetField(analyzer_in_filters, 0, "gamma_order");
    if (analyzer_in_gamma_order == NULL
	|| !mxIsNumeric(analyzer_in_gamma_order)
	|| mxIsComplex(analyzer_in_gamma_order)
	|| mxGetNumberOfElements(analyzer_in_gamma_order) != 1) {
      throw("Gfb_Analyzer_fprocess: 1st argument must be a Gfb_Analyzer\n"
	    "       structure, and must contain a structure field named\n"
	    "       `filters', containing the filter specifications.  Use\n"
	    "       the Gfb_Analyzer_new function to construct a valid\n"
	    "       Gfb_Analyzer structure");
    }
    size_t gamma_order = (size_t) mxGetScalar(analyzer_in_gamma_order);

    Filterbank_Data * fbd = Filterbank_Data_new(channels, gamma_order);
    m_filterbank_data = fbd; //< For proper deletion on exception.

    /*
     * extract filter coefficients, initial filter states, and normalization
     * factors from the matlab Gfb_Filter structures
     */
    for (size_t channel = 0; channel < channels; ++channel) {
      const mxArray * filter_type        = mxGetField(analyzer_in_filters,
						      channel,
						      "type");
      const mxArray * filter_coefficient = mxGetField(analyzer_in_filters,
						      channel,
						      "coefficient");
      const mxArray * filter_state       = mxGetField(analyzer_in_filters,
						      channel,
						      "state");
      const mxArray * filter_gamma_order = mxGetField(analyzer_in_filters,
						      channel,
						      "gamma_order");
      const mxArray * filter_norm        = mxGetField(analyzer_in_filters,
						      channel,
						      "normalization_factor");
      if (filter_type == NULL
	  || !mxIsChar(filter_type)
	  || mxGetM(filter_type) != 1
	  || mxGetN(filter_type) != strlen("Gfb_Filter")

	  || filter_coefficient == NULL
	  || !mxIsNumeric(filter_coefficient)
	  || !mxIsDouble(filter_coefficient)
	  || mxIsSparse(filter_coefficient)
	  || mxGetNumberOfElements(filter_coefficient) != 1

	  || filter_state == NULL
	  || !mxIsNumeric(filter_state)
	  || !mxIsDouble(filter_state)
	  || mxIsSparse(filter_state)
	  || mxGetM(filter_state) != 1
	  || size_t(mxGetN(filter_state)) != gamma_order
	
	  || filter_gamma_order == NULL
	  || !mxIsNumeric(filter_gamma_order)
	  || mxIsComplex(filter_gamma_order)
	  || mxGetNumberOfElements(filter_gamma_order) != 1
	  || mxGetScalar(filter_gamma_order) != gamma_order
	
	  || filter_norm == NULL
	  || !mxIsNumeric(filter_norm)
	  || mxIsComplex(filter_norm)
	  || mxGetNumberOfElements(filter_norm) != 1) {
	throw ("Gfb_Analyzer_fprocess: 1st argument must be a Gfb_Analyzer\n"
	       "       structure, and must contain a structure field named\n"
	       "       `filters', containing valid filter specifications.\n"
	       "       Use the Gfb_Analyzer_new function to construct a\n"
	       "       valid Gfb_Analyzer structure");
      }

      fbd->normalization_factors[channel] = mxGetScalar(filter_norm);

      fbd->real_filter_coefficients[channel] = *mxGetPr(filter_coefficient);
      /* permit real filter coefficients */
      fbd->imag_filter_coefficients[channel] =
	mxGetPi(filter_coefficient) ? *mxGetPi(filter_coefficient) : 0.0;

      FilterN * filter =
	new FilterN(fbd->gamma_order,
		    Complex(fbd->real_filter_coefficients[channel],
			    fbd->imag_filter_coefficients[channel]));

      size_t filter_stage;
      for (filter_stage = 0; filter_stage < gamma_order; ++filter_stage) {
	double state_r = mxGetPr(filter_state)[filter_stage];
	double state_i =
	  mxGetPi(filter_state) ? mxGetPi(filter_state)[filter_stage] : 0.0;
	
	fbd->real_filter_state[channel * gamma_order + filter_stage] = state_r;
	fbd->imag_filter_state[channel * gamma_order + filter_stage] = state_i;
	
	filter->stage(filter_stage).state = Complex(state_r, state_i);  
      }
      m_filter.push_back(filter);
    }
    return fbd;
  }

  void process(void) {
    Complex input;
    for (size_t sample = 0;
	 sample < m_filterbank_input->m_input_samples;
	 ++sample)
      for (size_t channel = 0; channel < m_channels; ++channel) {
	input.real = (m_filterbank_input->real(sample, channel) *
		      m_filterbank_data->normalization_factors[channel]);
	input.imag = (m_filterbank_input->imag(sample, channel) *
		      m_filterbank_data->normalization_factors[channel]);
	const Complex & output = m_filter[channel]->filter_complex(input);
	m_real_output[sample * m_channels + channel] = output.real;
	m_imag_output[sample * m_channels + channel] = output.imag;
      }

    /* copy analyzer struct if necessary */
    if (m_analyzer_out) {
      mxArray * filters = 0;
      *m_analyzer_out = mxDuplicateArray(m_analyzer_in);
      filters = mxGetField(*m_analyzer_out, 0, "filters");
      for (size_t channel = 0; channel < m_channels; ++channel) {
	mxArray * state_vector =
	  mxCreateDoubleMatrix(1, m_filterbank_data->gamma_order, mxCOMPLEX);
	unsigned filter_stage;
	for (filter_stage = 0;
	     filter_stage < m_filterbank_data->gamma_order;
	     ++filter_stage) {
	  mxGetPr(state_vector)[filter_stage] =
	    m_filter[channel]->stage(filter_stage).state.real;
	  mxGetPi(state_vector)[filter_stage] =
	    m_filter[channel]->stage(filter_stage).state.imag;
	}
	/* replace old filter state */
	mxArray * old_state_vector = mxGetField(filters, channel, "state");
	mxSetField(filters, channel, "state", state_vector);
	mxDestroyArray(old_state_vector);
      }
    }
  }
};
  



/**
 * The Matlab extension interface function. It checks for proper arguments,
 * extracts the needed data from the input arguments and creates the output
 * arguments.
 */
extern "C"
void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  try {
    Filterbank(nlhs, plhs, nrhs, prhs).process();
  }
  catch (const char * message) {
    mexErrMsgTxt(message);
  }
}
