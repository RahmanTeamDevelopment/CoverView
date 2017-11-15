"""
Maths utility functions.
"""


def min_or_nan(data):
    if len(data) == 0:
        return float('NaN')
    else:
        return min(data)


def max_or_nan(data):
    if len(data) == 0:
        return float('NaN')
    else:
        return max(data)
