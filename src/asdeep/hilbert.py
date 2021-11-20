"""
MIT License

Copyright (c) 2020 Princeton Laboratory for Intelligent Probabilistic Systems

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

# The code was copy-paste here to make the dependencies lighter
# Check the package: https://github.com/PrincetonLIPS/numpy-hilbert-curve

import numpy as np

def right_shift(binary, k=1, axis=-1):
  ''' Right shift an array of binary values.
  Parameters:
  -----------
   binary: An ndarray of binary values.
   k: The number of bits to shift. Default 1.
   axis: The axis along which to shift.  Default -1.
  Returns:
  --------
   Returns an ndarray with zero prepended and the ends truncated, along
   whatever axis was specified.
'''

  # If we're shifting the whole thing, just return zeros.
  if binary.shape[axis] <= k:
    return np.zeros_like(binary)

  # Determine the padding pattern.
  padding = [(0,0)] * len(binary.shape)
  padding[axis] = (k,0)

  # Determine the slicing pattern to eliminate just the last one.
  slicing = [slice(None)] * len(binary.shape)
  slicing[axis] = slice(None, -k)

  shifted = np.pad(binary[tuple(slicing)], padding,
                   mode='constant', constant_values=0)

  return shifted


def binary2gray(binary, axis=-1):
  ''' Convert an array of binary values into Gray codes.
  This uses the classic X ^ (X >> 1) trick to compute the Gray code.
  Parameters:
  -----------
   binary: An ndarray of binary values.
   axis: The axis along which to compute the gray code. Default=-1.
  Returns:
  --------
   Returns an ndarray of Gray codes.
  '''
  shifted = right_shift(binary, axis=axis)

  # Do the X ^ (X >> 1) trick.
  gray = np.logical_xor(binary, shifted)

  return gray


def decode(hilberts, num_dims, num_bits):
  ''' Decode an array of Hilbert integers into locations in a hypercube.

  This is a vectorized-ish version of the Hilbert curve implementation by John
  Skilling as described in:

  Skilling, J. (2004, April). Programming the Hilbert curve. In AIP Conference
    Proceedings (Vol. 707, No. 1, pp. 381-387). American Institute of Physics.

  Params:
  -------
   hilberts - An ndarray of Hilbert integers.  Must be an integer dtype and
              cannot have fewer bits than num_dims * num_bits.

   num_dims - The dimensionality of the hypercube. Integer.

   num_bits - The number of bits for each dimension. Integer.

  Returns:
  --------
   The output is an ndarray of unsigned integers with the same shape as hilberts
   but with an additional dimension of size num_dims.
  '''

  if num_dims*num_bits > 64:
    raise ValueError(
      '''
      num_dims=%d and num_bits=%d for %d bits total, which can't be encoded
      into a uint64.  Are you sure you need that many points on your Hilbert
      curve?
      ''' % (num_dims, num_bits)
    )

  # Handle the case where we got handed a naked integer.
  hilberts = np.atleast_1d(hilberts)

  # Keep around the shape for later.
  orig_shape = hilberts.shape

  # Treat each of the hilberts as a sequence of eight uint8.
  # This treats all of the inputs as uint64 and makes things uniform.
  hh_uint8 = np.reshape(hilberts.ravel().astype('>u8').view(np.uint8), (-1, 8))

  # Turn these lists of uints into lists of bits and then truncate to the size
  # we actually need for using Skilling's procedure.
  hh_bits = np.unpackbits(hh_uint8, axis=1)[:,-num_dims*num_bits:]

  # Take the sequence of bits and Gray-code it.
  gray = binary2gray(hh_bits)

  # There has got to be a better way to do this.
  # I could index them differently, but the eventual packbits likes it this way.
  gray = np.swapaxes(
    np.reshape(gray, (-1, num_bits, num_dims)),
    axis1=1, axis2=2,
  )

  # Iterate backwards through the bits.
  for bit in range(num_bits-1, -1, -1):

    # Iterate backwards through the dimensions.
    for dim in range(num_dims-1, -1, -1):

      # Identify which ones have this bit active.
      mask = gray[:,dim,bit]

      # Where this bit is on, invert the 0 dimension for lower bits.
      gray[:,0,bit+1:] = np.logical_xor(gray[:,0,bit+1:], mask[:,np.newaxis])

      # Where the bit is off, exchange the lower bits with the 0 dimension.
      to_flip = np.logical_and(
        np.logical_not(mask[:,np.newaxis]),
        np.logical_xor(gray[:,0,bit+1:], gray[:,dim,bit+1:])
      )
      gray[:,dim,bit+1:] = np.logical_xor(gray[:,dim,bit+1:], to_flip)
      gray[:,0,bit+1:] = np.logical_xor(gray[:,0,bit+1:], to_flip)

  # Pad back out to 64 bits.
  extra_dims = 64 - num_bits
  padded = np.pad(gray, ((0,0), (0,0), (extra_dims,0)),
                  mode='constant', constant_values=0)

  # Now chop these up into blocks of 8.
  locs_chopped = np.reshape(padded[:,:,::-1], (-1, num_dims, 8, 8))

  # Take those blocks and turn them unto uint8s.
  locs_uint8 = np.squeeze(np.packbits(locs_chopped, bitorder='little', axis=3))

  # Finally, treat these as uint64s.
  flat_locs = locs_uint8.view(np.uint64)

  # Return them in the expected shape.
  return np.reshape(flat_locs, (*orig_shape, num_dims))
