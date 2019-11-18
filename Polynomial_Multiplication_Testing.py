# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 07:29:33 2018

@author: Marcus
"""

from Polynomial_Multiplication import *

import pytest


def test_schoolbook_multiply():
    assert schoolbook_multiply([1,2],[1,2]) == [1,4,4]
schoolbook_multiply([2,3],[3,4])