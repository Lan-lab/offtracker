#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from importlib_metadata import distribution

dist = distribution('offtracker')
package_path = dist.locate_file('')
utility_dir = os.path.join(package_path, 'offtracker/utility')
os.chmod( os.path.join(utility_dir, 'bedGraphToBigWig'), 0o755)
print('offtracker is initialized.')