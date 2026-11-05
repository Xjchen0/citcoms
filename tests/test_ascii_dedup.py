import unittest

class TestAsciiDeduplication(unittest.TestCase):
    def parse_optional(self, opt_str):
        flags = {
            'connectivity': 0, 'stress': 0, 'pressure': 0, 'surf': 0,
            'botm': 0, 'geoid': 0, 'horiz_avg': 0, 'horiz_e2_avg': 0,
            'seismic': 0, 'coord_bin': 0, 'tracer': 0, 'comp_el': 0,
            'comp_nd': 0, 'strain_nd': 0, 'strain_el': 0, 'heating': 0,
            'viscy': 0, 'velo': 0, 'visc': 0, 'coord': 0
        }
        if not opt_str:
            return flags

        for token in opt_str.split(','):
            token = token.strip()
            if not token:
                continue
            if token == 'connectivity':
                flags['connectivity'] = 1
            elif token in ('velo', 'velocity'):
                flags['velo'] = 1
            elif token in ('visc', 'viscosity'):
                flags['visc'] = 1
            elif token in ('coord', 'coordinate'):
                flags['coord'] = 1
            elif token in ('3d', 'volume'):
                flags['velo'] = 1
                flags['visc'] = 1
                flags['coord'] = 1
            elif token == 'all':
                flags['velo'] = 1
                flags['visc'] = 1
                flags['coord'] = 1
                flags['surf'] = 1
                flags['botm'] = 1
                flags['horiz_avg'] = 1
                flags['stress'] = 1
                flags['pressure'] = 1
            elif token == 'horiz_avg':
                flags['horiz_avg'] = 1
            elif token == 'surf':
                flags['surf'] = 1
            elif token == 'botm':
                flags['botm'] = 1

        return flags

    def test_minimal_monitoring(self):
        # User's case: output_optional_ascii = "horiz_avg"
        flags = self.parse_optional("horiz_avg")
        self.assertEqual(flags['horiz_avg'], 1)
        self.assertEqual(flags['velo'], 0)
        self.assertEqual(flags['visc'], 0)
        self.assertEqual(flags['coord'], 0)

    def test_composite_ascii_h5(self):
        # Composite format "ascii,hdf5"
        # ASCII pass:
        opt_ascii = "horiz_avg"
        ascii_flags = self.parse_optional(opt_ascii)
        self.assertEqual(ascii_flags['velo'], 0, "ASCII must not output 3D velocity in composite mode")
        self.assertEqual(ascii_flags['visc'], 0, "ASCII must not output 3D viscosity in composite mode")
        self.assertEqual(ascii_flags['horiz_avg'], 1)

        # H5 pass:
        opt_h5 = "surf,botm,horiz_avg"
        h5_flags = self.parse_optional(opt_h5)
        self.assertEqual(h5_flags['surf'], 1)
        self.assertEqual(h5_flags['botm'], 1)
        self.assertEqual(h5_flags['horiz_avg'], 1)

    def test_explicit_3d_ascii(self):
        # User wants 3D volume in ASCII
        flags = self.parse_optional("velo,visc,horiz_avg")
        self.assertEqual(flags['velo'], 1)
        self.assertEqual(flags['visc'], 1)
        self.assertEqual(flags['horiz_avg'], 1)

        # Or using "3d" shorthand
        flags_3d = self.parse_optional("3d,horiz_avg")
        self.assertEqual(flags_3d['velo'], 1)
        self.assertEqual(flags_3d['visc'], 1)
        self.assertEqual(flags_3d['coord'], 1)
        self.assertEqual(flags_3d['horiz_avg'], 1)

if __name__ == '__main__':
    unittest.main()
