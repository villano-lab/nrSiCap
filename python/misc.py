#Miscellaneous Stuff

from math import log, floor

def human_format(number):
    units = ['', 'k', 'M', 'G', 'T', 'P']
    k = 1000.0
    magnitude = int(floor(log(number, k)))
    return '%.0f%s' % (number / k**magnitude, units[magnitude])