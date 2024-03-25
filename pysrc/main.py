#!/usr/bin/env python
import matplotlib.pyplot as plt
import Analysis

def main():
    data_dir = '../data/'

    comp_dir = data_dir + 'complexity/'
#    Analysis.AnalyseComplexity(comp_dir + 'complexity-clang.out',
#                               comp_dir + 'analysis/')
#    Analysis.AnalyseError(comp_dir + 'dump-clang/',
#                          comp_dir + 'error/')

#    p_comp_dir = data_dir + 'p-complexity/'
#    Analysis.AnalyseExpansionOrder(p_comp_dir + 'complexity.out',
#                                   p_comp_dir + 'analysis/')
#
#    Analysis.AnalyseExpansionError(p_comp_dir + 'dump/',
#                                   p_comp_dir + 'error/')
#    theta_comp_dir = data_dir + 'theta-complexity/'
#    Analysis.AnalyseTheta(theta_comp_dir + 'complexity.out',
#                          theta_comp_dir + 'analysis/')
#
#    Analysis.AnalyseThetaError(theta_comp_dir + 'dump/',
#                               theta_comp_dir + 'error/')

#    cold_dir = data_dir + 'cold/'
#    Analysis.AnalyseEvo(cold_dir, cold_dir + 'analysis/')

#    disk_dir = data_dir + 'disk/'
#    Analysis.AnalyseEvo(disk_dir, disk_dir + 'analysis/')

    two_galaxy_dir = data_dir + 'two-galaxies/'
    Analysis.AnalyseEvo(two_galaxy_dir, two_galaxy_dir + 'analysis/')


if __name__ == '__main__':
    main()
