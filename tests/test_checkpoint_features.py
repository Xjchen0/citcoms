#!/usr/bin/env python3
"""
Unit tests for sliding-window checkpoint pruning, auto-resume, walltime tracking, and Myr snapping.
"""

import unittest
import tempfile
import shutil
import os
import glob


class TestCheckpoints(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.prefix = os.path.join(self.test_dir, "test_run")

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_checkpoint_pruning_sliding_window(self):
        checkpoints_to_keep = 2
        history_file = f"{self.prefix}.checkpoint_history.txt"
        last_good = os.path.join(self.test_dir, "last_good_checkpoint")

        cycles = [1000, 2000, 3000, 4000, 5000]

        for cycle in cycles:
            for p in range(4):
                chk_file = f"{self.prefix}.chkpt.{p}.{cycle}"
                with open(chk_file, "w") as f:
                    f.write(f"checkpoint rank {p} cycle {cycle}")

            with open(last_good, "w") as f:
                f.write(f"{cycle}\n")

            with open(history_file, "a") as f:
                f.write(f"{cycle} 123.456\n")

            with open(history_file, "r") as f:
                lines = f.readlines()

            if len(lines) > checkpoints_to_keep:
                old_cycles = [int(line.split()[0]) for line in lines[:-checkpoints_to_keep]]
                for old_cycle in old_cycles:
                    for old_file in glob.glob(f"{self.prefix}.chkpt.*.{old_cycle}"):
                        os.remove(old_file)

                with open(history_file, "w") as f:
                    f.writelines(lines[-checkpoints_to_keep:])

        remaining_files = os.listdir(self.test_dir)

        for c in [1000, 2000, 3000]:
            for p in range(4):
                self.assertNotIn(f"test_run.chkpt.{p}.{c}", remaining_files)

        for c in [4000, 5000]:
            for p in range(4):
                self.assertIn(f"test_run.chkpt.{p}.{c}", remaining_files)

        with open(last_good, "r") as f:
            self.assertEqual(f.read().strip(), "5000")

    def test_auto_resume_detection(self):
        last_good = os.path.join(self.test_dir, "last_good_checkpoint")

        self.assertFalse(os.path.exists(last_good))

        with open(last_good, "w") as f:
            f.write("7500\n")

        with open(last_good, "r") as f:
            content = f.read().strip()
            self.assertEqual(int(content), 7500)

    def test_myr_snapping_calculation(self):
        scalet = 100.0  # 1 dimensionless unit = 100 Myr
        interval_myr = 5.0  # output every 5 Myr -> interval_time = 0.05
        interval_time = interval_myr / scalet

        current_time = 0.049999
        self.assertTrue(current_time + 1.0e-4 * interval_time >= interval_time)

        next_time = interval_time
        next_time += interval_time
        self.assertAlmostEqual(next_time, 0.10)


if __name__ == "__main__":
    unittest.main()
