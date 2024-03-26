#!/usr/bin/env python
import matplotlib.pyplot as plt
import Analysis

def main():
    data_dir = '../data/'

    n_comp_dir = data_dir + 'n-complexity/'
#    Analysis.AnalyseN(n_comp_dir + 'complexity.out',
#                      n_comp_dir + 'analysis/')
#    Analysis.AnalyseParamError(n_comp_dir + 'dump/',
#                               n_comp_dir + 'error/',
#                               paramName = 'n',
#                               xlabel = 'number of particles $n$')
#
#    p_comp_dir = data_dir + 'p-complexity/'
#    Analysis.AnalyseP(p_comp_dir + 'complexity.out',
#                                   p_comp_dir + 'analysis/')
#
#    Analysis.AnalyseParamError(p_comp_dir + 'dump/',
#                               p_comp_dir + 'error/',
#                               paramName = 'p',
#                               xlabel = 'order of expansion $p$')
#    theta_comp_dir = data_dir + 'theta-complexity/'
#    Analysis.AnalyseTheta(theta_comp_dir + 'complexity.out',
#                          theta_comp_dir + 'analysis/')
#
#    Analysis.AnalyseParamError(theta_comp_dir + 'dump/',
#                               theta_comp_dir + 'error/',
#                               paramName = 'theta',
#                               xlabel = 'opening angle $\\theta$')
#
#    cold_dir = data_dir + 'cold/'
#    Analysis.AnalyseEvo(cold_dir + 'dump/',
#                        cold_dir + 'analysis/')

    disk_dir = data_dir + 'disk/'
    Analysis.AnalyseEvo(disk_dir + 'dump/',
                        disk_dir + 'analysis/')
#
#    galaxy_dir = data_dir + 'galaxy/'
#    Analysis.AnalyseEvo(galaxy_dir + 'dump/',
#                        galaxy_dir + 'analysis/')
#
#    two_galaxy_dir = data_dir + 'two-galaxies/'
#    Analysis.AnalyseEvo(two_galaxy_dir + 'dump/',
#                        two_galaxy_dir + 'analysis/')
#

if __name__ == '__main__':
    main()
