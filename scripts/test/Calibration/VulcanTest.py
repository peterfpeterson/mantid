# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +

import unittest
from mantid.simpleapi import DeleteWorkspace, LoadEmptyInstrument
import numpy as np
#from mantid.simpleapi import CompareWorkspaces, DeleteWorkspaces, LoadNexusProcessed, UnGroupWorkspace
from Calibration.vulcan import make_group_workspace


class TestCreateGroupingWorkspace(unittest.TestCase):
    PIXELS_PER_TUBE = 512
    TUBES_PER_PANEL = 8
    PIXELS_PER_PANEL = PIXELS_PER_TUBE * TUBES_PER_PANEL
    BANK1_PANELS = 20
    BANK2_PANELS = 20
    BANK5_PANELS = 9
    TOTAL_TUBES = TUBES_PER_PANEL * (BANK1_PANELS + BANK2_PANELS + BANK5_PANELS)
    TOTAL_PIXELS = PIXELS_PER_PANEL * (BANK1_PANELS + BANK2_PANELS + BANK5_PANELS)
    instr = ""  # empty means unset

    @classmethod
    def setUpClass(cls) -> None:
        FILENAME = 'VULCAN_Definition.xml'
        instr = LoadEmptyInstrument(Filename=FILENAME, OutputWorkspace='VULCAN')
        assert instr, f'Empty instrument from {FILENAME}'
        cls.instr = str(instr)

    @classmethod
    def tearDownClass(cls) -> None:
        if cls.instr:
            DeleteWorkspace(cls.instr)
            cls.instr = ""

    def __create_grouping_plan(self, groups_per_tube: int):
        bank1_start = 0
        bank2_start = self.BANK1_PANELS * self.PIXELS_PER_PANEL
        bank5_start = (self.BANK1_PANELS + self.BANK2_PANELS) * self.PIXELS_PER_PANEL
        bank5_stop = self.TOTAL_PIXELS
        return [(bank1_start, groups_per_tube, bank2_start),
                (bank2_start, groups_per_tube, bank5_start),
                (bank5_start, groups_per_tube, bank5_stop)]

    def __create_unique_groups(self, groups_per_tube: int):
        # this is a (somewhat) silly way to designate that a tube strategy is wanted
        if groups_per_tube > 0:
            tube_strategy = self.__create_grouping_plan(groups_per_tube)
        else:
            tube_strategy = []

        grouping = make_group_workspace(self.instr, 'make_simple', group_detectors_by='bank', grouping_plan=tube_strategy)
        groups = np.unique(grouping.extractY())
        grouping.delete()  # workspace no longer needed

        return groups

    def __verify_groups(self, groups, num_groups_exp):
        assert len(groups) == num_groups_exp, 'Found {} != {} groups: {}'.format(num_groups_exp, len(groups), groups)
        assert np.min(groups) == 1, 'Expected group minimum = 1, found {}'.format(int(np.min(groups)))
        assert np.max(groups) == num_groups_exp, 'Expected group maximum = {}, found {}'.format(num_groups_exp,
                                                                                                int(np.max(groups)))

    def test_make_simple(self):
        groups = self.__create_unique_groups(0)
        self.__verify_groups(groups, 3)

    def test_make_full_tube(self):
        groups = self.__create_unique_groups(512)
        self.__verify_groups(groups, self.TOTAL_TUBES)  # one group per tube

    def test_make_half_tube(self):
        groups = self.__create_unique_groups(256)
        self.__verify_groups(groups, 2 * self.TOTAL_TUBES)  # two group per tube


if __name__ == '__main__':
    unittest.main()
