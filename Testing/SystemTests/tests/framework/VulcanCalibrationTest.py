# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +

import unittest
import os
import systemtesting
import tempfile
from Calibration.vulcan import load_and_crop, cross_correlate_calibrate, align_vulcan_data, peak_position_calibrate
#from vulcan.calibration.calibrate_vulcan_x import align_vulcan_data, cross_correlate_calibrate, peak_position_calibrate


class VulcanCalibrationTest(systemtesting.MantidSystemTest): # unittest.TestCase):
    diamond_files = ["VULCAN_192226.nxs.h5"]

    def requiredFiles(self):
        return self.diamond_files

    def runTest(self):
        test_output_dir = tempfile.gettempdir()
        calibration_rst = os.path.join(test_output_dir, "VULCAN_Calibration_Hybrid.h5")
        tube_grouping_plan = None

        # load the data
        diamond_ws = load_and_crop(runnumbers=self.diamond_files, bin_step=-0.0003, bad_pulse_threashold=0, convert_to_dspace=True,
                                   user_idf='VULCAN')

        # perform calibration
        # step_1:
        cc_calib_file, _ = cross_correlate_calibrate(diamond_ws, cross_correlate_param_dict=None,
                                                     output_dir=test_output_dir)
        print('CC CALIBRATION:', cc_calib_file)

        # step_2:
        cc_focus_ws_name = align_vulcan_data(
            diamond_runs=[diamond_ws],  # TODO second round should go from file
            diff_cal_file_name=cc_calib_file,
            output_dir=test_output_dir,
            tube_grouping_plan=tube_grouping_plan)

        # step_3:
        tube_grouping_plan = [(0, None, 81920), (81920, None, 81920 * 2), (81920 * 2, None, 200704)]
        peak_position_calibrate(cc_focus_ws_name, tube_grouping_plan, cc_calib_file, calibration_rst,
                                test_output_dir)

        # verify calibration file is created
        # NOTE: this is step 1, we will add more reliable asseration once
        #       we have a reference calibration results
        assert os.path.exists(calibration_rst)

    #def validateMethod(self):
    #    return None # it running is all that we need

    #def validate(self):
    #    pass


'''OLD prototype code
dia_runs = [192226, 192227] #, 192228, 192229, 192230]
diamond_ws = load_and_crop(runnumbers=dia_runs, bin_step=-0.0003, bad_pulse_threashold=0, convert_to_dspace=True,
                           user_idf='VULCAN')
#cross_correlate_calibrate called in calibrate_vulcan_x.pyL298
tube_cc_plan = cross_correlation_in_tubes()
cc_calib_file, diamond_ws_name = cross_correlate_calibrate(diamond_ws, cross_correlate_param_dict = tube_cc_plan)
print('file:', cc_calib_file)
#create grouping for sub-sections of the instrument
tube_grouping_plan = [(0, 512, 81920), (81920, 1024, 81920 * 2), (81920 * 2, 256, 200704)]  # TODO
#AlignAndFocusPowder using that information
cc_focus_ws_name = align_vulcan_data(diamond_runs=diamond_ws_name,
                                                     diff_cal_file_name=cc_calib_file,
                                                     output_dir='.', #output_dir,
                                                     tube_grouping_plan=tube_grouping_plan,
                                                     user_idf='VULCAN')

#peak_position_calibrate in calibrate_vulcan_x.pyL319 (does this combine the two parts of the DIFC calculation?)
final_calib_file = 'VULCAN_Calibration_Hybrid.h5'
peak_position_calibrate(cc_focus_ws_name, tube_grouping_plan, cc_calib_file, final_calib_file, output_dir='.')

#create the actual grouping workspace
#group_ws = 'VULCAN_group_banks'
#group_ws = CreateGroupingWorkspace(InputWorkspace=cc_focus_ws_name,
#                                   GroupDetectorsBy='bank',
#                                   OutputWorkspace=group_ws)
#SaveDIFC
'''
if __name__ == "__main__":
    unittest.main()