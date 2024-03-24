import Analysis

def main():
    Analysis.AnalyseComplexity('../data/complexity/complexity-clang.out',
                               '../data/complexity/analysis/')
    Analysis.AnalyseError('../data/complexity/dump-clang/',
                          '../data/complexity/analysis/error/')

if __name__ == '__main__':
    main()
