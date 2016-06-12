import os
import sys
import _04_Pmf
import _05_myplot
import _13_Cdf


def ReadTimeLength():
    filename = os.path.join('.', 'babyboom.csv')
    fp = open(filename)
    result = []

    for i, line in enumerate(fp):
        line = line.strip()
        result.append(int(line.split(',')[3]))

    fp.close()
    return result


def CalcTimeInterval():
    timeLength = ReadTimeLength()
    print('timeLength\'s length=', len(timeLength))
    print('range(len(timeLength))-->', range(len(timeLength)))
    timeInterval = [timeLength[i + 1] - timeLength[i] for i in range(len(timeLength) - 1)]
    return timeInterval


def CdfTimeInterval():
    timeInterval = CalcTimeInterval()
    pmf = _04_Pmf._make_pmf_from_list(timeInterval, 'baby birth interval')
    _05_myplot._clf()
    _05_myplot._hist(pmf)
    _05_myplot._save(root='baby_birth_interval_pmf',
                     title='PMF of baby birth interval',
                     xlabel='interval(minutes)',
                     ylabel='probability')

    cdf = _13_Cdf.MakeCdfFromPmf(pmf)

    _05_myplot._clf()
    _05_myplot._cdf(cdf)
    _05_myplot._save(root='baby_birth_interval_cdf',
                     title='CDF of baby birth interval',
                     xlabel='interval(minutes)',
                     ylabel='cumulative probability')


def LogCCdfTimeIntervel():
    timeInterval = CalcTimeInterval()
    pmf = _04_Pmf._make_pmf_from_list(timeInterval, 'baby birth interval')
    cdf = _13_Cdf.MakeCdfFromPmf(pmf, 'baby birth interval')
    _05_myplot._clf()
    _05_myplot._cdf(cdf, complement=True, xscale='linear', yscale='log')
    _05_myplot._save(root='baby_birth_interval_logccdf',
                     title='LogCCDF of baby birth interval',
                     xlabel='interval(minutes)',
                     ylabel='LogCCdf')


# print CalcTimeInterval()
# print ReadTimeLength()
# CdfTimeInterval()
# print sorted(CalcTimeInterval())
LogCCdfTimeIntervel()
