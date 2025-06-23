#!/usr/bin/env python3

import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
sys.path.append(parent_dir)

from callLlm import call_llm
from callSplice import call_splice

spliceai_data = call_splice("chr8:140300616:T:G")
print("SpliceAI Data:")
print(spliceai_data)