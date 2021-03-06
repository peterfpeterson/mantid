.. algorithm::

.. summary::

.. relatedalgorithms::

.. properties::

Description
-----------

Reads a McStas Nexus file into a Mantid WorkspaceGroup with the name of 
the algorithm property OutputWorkspace. Data generated by McStas monitor components are 
stored in workspaces of type Workspace2D and/or EventWorkspace. The name of a
workspace equals that of the mcstas component name + '_' + name of the OutputWorkspace.
In addition an EventWorkspace with the name 'EventData' + '_' + name of the OutputWorkspace
is created which contains the sum of all event datasets in the Nexus file.
Note if OutputOnlySummedEventWorkspace=True only this EventWorkspace is returned
by the algorithm.

For information about how to create McStas outputs that can 
readily be read by this loader, see `here <https://github.com/McStasMcXtrace/McCode/wiki/McStas-and-Mantid>`_.
For more information about McStas, and combined McStas and Mantid analyses see references below.

The ErrorBarsSetTo1 property applies to event data, but not to histogram data.

LoadMcStas replaces LoadMcStasEventNexus. LoadMcStas can be used for 
reading McStas 2.1 histogram and event data. 
:ref:`algm-LoadMcStasNexus` can be used 
for reading McStas 2.0 histogram data. 

Information about the required structure of the input Nexus file
################################################################

The input file must have an 
``/entry1/simulation/name`` item whose value is ``"mccode"``.

The output workspace will contain one workspace for each group of
class ``NXdata`` in the input NeXus file, which is not of name ``"content_nxs"`` 
found in a group of class ``NXDetector`` of name ``"Data"``.
The name of the workspace is the same as the name of the group, 
but with with the name of the outputworkspace, as specified in the algorithm properties,
added to the end.

+----------------------------------+------------------------------------------+---------------------------------------+
| Description of Data              | Found in Nexus file (within 'run')       | Placed in Workspace (Workspace2D)     |
|                                  |                                          | or output                             |
+==================================+==========================================+=======================================+
| All data                         | Group of class ``NXDetector``            | See below                             |
|                                  | of name ``"data"``                       |                                       |
+----------------------------------+------------------------------------------+---------------------------------------+
| Generic group including either   | each group of class ``NXdata`` in        | one workspace each                    |
| event or histrogram data         | ``"data"``, henceforth referred to as    |                                       |
|                                  | [DATA]                                   |                                       |
+----------------------------------+------------------------------------------+---------------------------------------+
| Event data                       | item in a [DATA] with a ``long_name``    | event data                            |
|                                  | attribute containing ``"Neutron_ID"``    |                                       |
|                                  | and name ``"events"``                    |                                       |
+----------------------------------+------------------------------------------+---------------------------------------+
| Histrogram data                  | items in a [DATA] with a ``long_name``   | histogram data                        |
|                                  | attribute **not** containing             |                                       |
|                                  | ``"Neutron_ID"``                         |                                       |
+----------------------------------+------------------------------------------+---------------------------------------+
| Instrument                       | ``/instrument``                          | loaded into workspace, only if        |
|                                  |                                          | events are loaded                     | 
+----------------------------------+------------------------------------------+---------------------------------------+
| Instrument definition            | ``/instrument/instrument_xml/data``      | loaded into workspace, only if        |
|                                  | is needed for events to be loaded        | events are loaded                     | 
+----------------------------------+------------------------------------------+---------------------------------------+

The event data of the McStas file occurs in a NeXus table with six columns:

1. Weight
2. X coordinate
3. Y coordinate
4. Number of previous events
5. Detector ID
6. Time


References
##########

For more information about McStas and its general usage for simulating neutron 
scattering instruments and experiments visit the `McStas homepage <http://www.mcstas.org>`_ .

For examples of how combined McStas and Mantid analyses can help 
instrument simulation and data treatment/analysis tasks see Nielsen., T.R. et al., McStas
and Mantid integration, Journal of Neutron Research, vol. 18, no. 2-3, pp. 61-77, 2015
`DOI: 10.3233/JNR-160026 <http://dx.doi.org/10.3233/JNR-160026>`_ 
[`arXiv <http://arxiv.org/abs/1607.02498>`_].

Usage
-----

.. include:: ../usagedata-note.txt

**Example - Load McStas data containing both event and histogram data:**

.. testcode:: ExLoadMcStas

   # Load the data into tuple
   ws = LoadMcStas('mcstas_event_hist.h5')

   # workspace group is first entry in tuple
   group = mtd['ws']
   print("Number of entries in group: {}".format(group.getNumberOfEntries()))

   eventData = mtd['EventData_ws']
   print("Number of histograms in event data: {}".format(eventData.getNumberHistograms()))
   print("Name of event data: {}".format(eventData.getName()))

   someHistogramData = mtd['Edet.dat_ws']
   print("Number of histograms in hist data: {}".format(someHistogramData.getNumberHistograms()))
   print("Name of hist data: {}".format(someHistogramData.getName()))


Output:

.. testoutput:: ExLoadMcStas

   Number of entries in group: 5
   Number of histograms in event data: 8192
   Name of event data: EventData_ws
   Number of histograms in hist data: 1
   Name of hist data: Edet.dat_ws

**Example - Comparing event data entries in a McStas Nexus file:**

The mccode_multiple_scattering.h5 McStas Nexus file contains two event data entries:
named single_list_p_x_y_n_id_t and multi_list_p_x_y_n_id_t, one
from each of two detector banks of the instrument simulated. Setting 
OutputOnlySummedEventWorkspace=False these are loaded
individually into separate workspaces. In addition, this algorithm returns the
workspace EventData_ws, which contains the sum of all event data entries in the
McStas Nexus file. The example below performs a test to show that the summation
of the workspaces has been executed correctly.

.. testcode:: CheckEqualScattering

    # Load the data into tuple
    ws = LoadMcStas('mccode_multiple_scattering.h5', OutputOnlySummedEventWorkspace=False)

    # Calculate total of all event data entries
    all_scattering_event_ws = mtd['EventData_ws']
    total_all = 0
    for i in range(all_scattering_event_ws.getNumberHistograms()):
      total_all += all_scattering_event_ws.readY(i)[0]
    print("The sum of all scattering spectra: {0:.6e}".format(total_all))

    # Calculate total scattering from the single event bank
    single_scatter_event_ws = mtd['single_list_p_x_y_n_id_t_ws']
    total_single = 0
    for i in range(single_scatter_event_ws.getNumberHistograms()):
      total_single += single_scatter_event_ws.readY(i)[0]
    print("The sum of all single scattering spectra: {0:.6e}".format(total_single))

    # Calculate total scattering from the 'k02' detector bank
    multiple_scatter_event_ws = mtd['multi_list_p_x_y_n_id_t_ws']
    total_multiple = 0
    for i in range(multiple_scatter_event_ws.getNumberHistograms()):
      total_multiple += multiple_scatter_event_ws.readY(i)[0]
    print("The sum of all multiple scattering spectra: {0:.6e}".format(total_multiple))

    # Check equality
    sum_of_scattering = total_multiple + total_single
    # This is equal to the sum of all scattering spectra
    print("Sum of single and multiple scattering workspaces: {0:.6e}".format(total_single + total_multiple))

Output:

.. testoutput:: CheckEqualScattering

   The sum of all scattering spectra: 2.038678e-11
   The sum of all single scattering spectra: 1.907862e-11
   The sum of all multiple scattering spectra: 1.308161e-12
   Sum of single and multiple scattering workspaces: 2.038678e-11

.. categories::

.. sourcelink::
