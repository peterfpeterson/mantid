.. _Powder Diffractions Sample Corrections:

Powder Diffractions Sample Corrections
======================================

.. contents::

The work to be done
-------------------

It will start with some equations from Appendix B

describe the full instrument -> partially grouped / low resolution
instrument (need a name) -> final spectra (6 for NOMAD Rietveld, 1 for
everything else)

show how the equations become simpler when you ignore various sample
corrections, and show how that compares (mathematically) to what
SNSPowderReduction does currently.

That is something that can be verified by Joerg / Thomas / you /
???. Then the document will compare the various absorption
corrections.

The reason for all the equations at the beginning is to 1. not forget
what I finally understood today 2. give context to how absorption fits
into absolute normalization. Then the big bit of thinking left to do
is sorting out how becomes a mantid workflow algorithm that supports
verifying various corrections (for auditing), re-reducing data with
tweaks to parameters, and event filtering with summing.

Full correction equations
-------------------------

To generate a reciprocal space spectrum for use in total scattering
(i.e. to Fourier transform) one needs to get to the differential cross
section (DCS) and scale it with the coherent scattering length (note
that that these subscripts are described :ref:`here <Sample
Corrections>`).

.. math::

   S(Q) = \frac{1}{\langle b_{coh} \rangle^2} \frac {d \sigma}{d\Omega}
          - \frac{\langle b^2_{tot}\rangle - \langle b_{coh} \rangle^2}{\langle b_{coh} \rangle^2}

The second part of this is called the normalized Laue term. The DCS is
defined as

.. math::
   :label: the_big_one

   \frac{d\sigma}{d\Omega} =   \frac{1}{A_{s,sca} N_s}   \bigg [ \frac{1}{\Phi} (I^E_s - I^E_B) - M_{sca}
                               - \frac{A_{c,sca}}{A_{c,ca}} \Big [\frac{1}{\Phi} (I^E_c - I^E_B) - M_{ca}
                               - \frac{A_{a,ca}}{A_{a,a}} [\frac{1}{\Phi} (I^E_a - I^E_B) - M_{a}] \Big ]
                               - \frac{A_{a,sca}}{A_{a,a}} \Big [ \frac{1}{\Phi} (I^E_a - I^E_B) - M_{a} \Big ] \bigg ]  - P_s^{ie}

In this equation, :math:`M` represents multiple scattering, :math:`A`
is the absorption, and :math:`P` is the inelastic scattering. If the
empty sample environment is not measured, then those terms drop
out. The normalization term, :math:`\Phi`, is derived from the
vanadium being processed in a way similar to Eq. :eq:`the_big_one`.
making the assumption that differential cross section for vanadium
:math:`\frac{d\sigma}{d\omega} \approx \frac{\sigma_{tot}}{4\pi} =
b_{tot}^2`.

.. math::
   :label: the_big_one_simplified_for_vanadium

   \Phi = \frac{ I^E_n - I^E_B } {[A_{n,n} N_n [ \langle b_v^2 \rangle + P_n^{ie} ] + M_{n}]  }

As you can see, this relies heavily on the algorithms described in :ref:`Sample Corrections`.

Various simplifications
-----------------------

Dropping the sample environment term, and making the conversion of
:math:`\frac{A_{c,sc}}{A_{c,c}A_{s,sc}} = \frac{1}{A_{c,c}}` and
:math:`\frac{1}{A_{s,sc}} = \frac{1}{A_{s,s}}` for when partial
absorption calculations weren't made, Eq. :eq:`the_big_one` becomes

.. math::

   \frac{d\sigma}{d\Omega} =   \frac{1}{A_{s,s} N_s}   \bigg [ \frac{1}{\Phi} (I^E_s - I^E_B) - M_{sc}
                               - \frac{A_{s,s}}{A_{c,c}} \Big [\frac{1}{\Phi} (I^E_c - I^E_B) - M_{c} \Big ] \bigg ]  - P_s^{ie}

with the furnace correction

.. math::

   \frac{d\sigma}{d\Omega} =   \frac{1}{A_{s,s} N_s} \bigg [ \frac{1}{\Phi} (I^E_s - I^E_B) - M_{sca}
                               - \frac{A_{s,s}}{A_{c,c}} \Big [\frac{1}{\Phi} (I^E_c - I^E_B) - M_{ca}
                               - \frac{A_{c,c}}{A_{a,a}} [\frac{1}{\Phi} (I^E_a - I^E_B) - M_{a}] \Big ]
                               - \frac{A_{s,s}}{A_{a,a}} \Big [ \frac{1}{\Phi} (I^E_a - I^E_B) - M_{a} \Big ] \bigg ]  - P_s^{ie}


In the case of not having multiple scattering, or when the decision
was made that it is not big enough to warrent correcting, the
normalization factor can be brought out front and the absorption can
be applied pixel-by-pixel before grouping into super-pixels (read
discription below) and having the inelastic scattering correction
applied.

.. math::

   \frac{d\sigma}{d\Omega} = \frac{1}{\Phi N_s} \frac{1}{A_{s,sca}}   \bigg [ (I^E_s - I^E_B)
                               - \frac{A_{c,sca}}{A_{c,ca}} \Big [(I^E_c - I^E_B)
                               - \frac{A_{a,ca}}{A_{a,a}} [ (I^E_a - I^E_B)] \Big ]
                               - \frac{A_{a,sca}}{A_{a,a}} \Big [(I^E_a - I^E_B) \Big ] \bigg ]  - P_s^{ie}

The even bigger simplification of no multiple scattering or empty
sample environment which contains the kernel of the equation in
:ref:`ApplyPaalmanPingsCorrection <algm-ApplyPaalmanPingsCorrection>`.

.. math::

   \frac{d\sigma}{d\Omega} = \frac{1}{\Phi N_s} \frac{1}{A_{s,sc}}   \bigg [ (I^E_s - I^E_B)
                               - \frac{A_{c,sc}}{A_{c,c}} \Big [(I^E_c - I^E_B)
                               \Big ] \bigg ]  - P_s^{ie}

Finally, ignoring multiple scattering, inelastic scattering, and
assuming that only the full absorption corrections are calculated
yields

.. math::

   \frac{d\sigma}{d\Omega} =   \frac{1}{\Phi N_s} \bigg [ \frac{1}{A_{s,s}} (I^E_s - I^E_B)
                                                        - \frac{1}{A_{c,c}} (I^E_c - I^E_B) \bigg ]

The careful reader will notice that this is similar to what is
currently implemented in :ref:`SNSPowderReduction
<algm-SNSPowderReduction>` minus the factor of :math:`N_s` and with
the normalization created slightly differently. That code is in the
process of being upgraded to the equation (lots of detail missing and
:math:`I(Q)` being an overloaded term in this case)

.. math::

   I(Q) = \bigg [ \frac{1}{A_{s,s}} (I^E_s - I^E_B) - \frac{1}{A_{c,c}} (I^E_c - I^E_B) \bigg ]
          / \bigg [ \frac{1}{A_{v,v}} (I^E_v - I^E_B) \bigg ]

Using this fancy equation in practice
-------------------------------------

While all of these correction can, in principle, be applied at the
pixel-by-pixel level, characterization measurements (e.g. vanadium,
empty container, empty instrument) are often lacking sufficient
statistics to allow doing so without the uncertainties from the
characterizations overwhelming the uncertainties in the final
:math:`S(Q)`. This must be balanced with the various corrections
losing accuracy by summing detector pixels together. For this reason,
we introduce the concept of "partially grouped/super-pixels/low
resolution instrument/need a name." With this concept in hand, we can
enumerate how one can implement Eq. :eq:`the_big_one`.

1. :ref:`AlignAndFocusPowderFromFiles
   <algm-AlignAndFocusPowderFromFiles>` to load the instrument, do
   various bits of summing and filtering and reduce the instrument to
   super-pixels. This should be done to all of the data involved in
   the overall reduction: sample, empty container, empty sample
   environment, vanadium rod, empty instrument, ... If one or more
   were not measured then various bits of this get simplier.
2. Subtract the empty instrument, :math:`I^E_B`, from all of the
   signals: sample (:math:`I^E_S`), empty container (:math:`I^E_C`),
   and empty sample environment :math:`I^E_a`).
3. Scale the super-pixels by the incident spectrum, :math:`\Phi`. *cache?*
4. Calculate the :ref:`multiple scattering <Multiple Scattering
   Corrections>` for the super-pixels and subtract it from the various
   signals listed in step 2.
5. Subtract the multiple scattering terms. *cache?*
6. Calculate the :ref:`absorption terms <Absorption
   Corrections>`. Since Eq. :eq:`the_big_one` uses the partials, it
   will need to be simplified to use full absorption terms. The
   relation of the two equations is seen in
   :ref:`ApplyPaalmanPingsCorrection
   <algm-ApplyPaalmanPingsCorrection>`. **The big equation needs a
   version that doesn't rely on partials since they are, largely,
   unavailable.**
7. Apply the various absorption terms to the square bracket parts of
   Eq. :eq:`the_big_one`. *cache?*
8. Do the top level subtraction
9. Apply the final absorption workspace and divide by number of atoms in the beam.
10. Calculate the inelastic scattering, the Plazcek correction will do
    **There isn't an algorithm yet**
10. Subtract the inelastic scattering
11. Divide the DCS by :math:`\langle b_{coh} \rangle^2` and subtract off the
    normalized Laue term.
12. :ref:`MatchSpectra <algm-MatchSpectra>` the super-pixels. This
    would probably be easier to use with :math:`S(Q)-1` **shouldn't be
    necessary**
13. :ref:`SumSpectra <algm-SumSpectra>` with ``WeightedSum=True,
    MultiplyBySpectra=False`` to get a single output workspace that is
    :math:`S(Q)` **missing a good way to create 6 spectra for Rietveld
    from the super-pixels**
14. Profit!!!!!!!!!!

Processing the vanadium would/does go through a similar process, with
the additions of :ref:`removing the Bragg peaks
<algm-StripVanadiumPeaks>` and smoothing the data. Creating a workflow
algorithm for processing vanadium would be useful. It is not uncommon
to smooth the container and emtpy instrument measurements as well.

The biggest thing to notice from the previous section on
simplifications of the equations, is that droping the multiple
scattering term allows for applying the absorption correction on a
pixel-by-pixel basis. This suggests the utility of determining whether
the correction is necessary. There should also be some effort to
determine various quantities that could be cached for reuse in other
reductions. Currently it is just the result of
:ref:`AlignAndFocusPowderFromFiles
<algm-AlignAndFocusPowderFromFiles>`, but the fully processed
normalization, :math:`\Phi` would greatly improve performance.
