import Analysis

def main():
    Analysis.AnalyseComplexity('../data/complexity/complexity-clang.out',
                               '../data/complexity/analysis/')
#    Analysis.AnalyseError('../data/complexity/dump-clang/',
#                          '../data/complexity/analysis/error/')
    Analysis.AnalyseExpansionOrder('../data/complexity/p-comp.out',
                                   '../data/complexity/analysis/p/')


if __name__ == '__main__':
    main()
