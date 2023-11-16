import numpy as np
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--input", dest="input", default='Bias_mu.txt', type="str", help="input file")
parser.add_option("-c", "--channel", dest="channel", default='mu', type="str", help="ele/mu?")
(options, args) = parser.parse_args()

def main():
    year=['16','17','18']
    mass=['M1','M2','M3','M4','M5','M6','M7','M8','M9','M10','M15','M20','M25','M30']
    bias_test={}
    bias_name={}
    bias_true={}
    
    f = open(options.input)
    f_lines = f.readlines()

    for m in range(len(mass)):
        bias_test[mass[m]]=[]
        bias_name[mass[m]]=[]
        bias_true[mass[m]]=[]
        exp = f_lines[3*m].rstrip('\n').rstrip().split('\t')
        exp_name = f_lines[3*m+1].rstrip('\n').split('\t')
        for i in range(len(exp)):
            bias_test[mass[m]].append(float(exp[i]))
            bias_name[mass[m]].append(exp_name[i])
        bias_true[mass[m]].append(round(float(f_lines[3*m+2].rstrip('\n')),2))

    '''
    #
    bias_test_M1=[-0.06018,-0.01283,0.0372,-0.01167]#Bias for each test function
    bias_name_M1=['1st Order Bernstein', '1st Order Exponential', '1st Order Laurent', '1st Order Power Low']
    bias_true_M1=[12.74]

    bias_test_M5=[-0.06026,-0.0469,0.014,-0.01128,-0.01944]#Bias for each test function
    bias_name_M5=['1st Order Bernstein', '2nd Order Bernstein', '1st Order Exponential', '1st Order Laurent', '1st Order Power Low']
    bias_true_M5=[1.95]

    bias_test_M15=[-0.05261,-0.0217,-0.02854,-0.01751,-0.004545]#Bias for each test function
    bias_name_M15=['1st Order Bernstein', '2nd Order Bernstein', '1st Order Exponential', '1st Order Laurent', '1st Order Power Low']
    bias_true_M15=[4.11]

    bias_test_M30=[0.004876,-0.02384,-0.03709,0.005416]#Bias for each test function
    bias_name_M30=['1st Order Bernstein', '1st Order Exponential', '1st Order Laurent', '1st Order Power Low']
    bias_true_M30=[2.28]

    bias_test_all=[bias_test_M1, bias_test_M5, bias_test_M15, bias_test_M30]
    bias_name_all=[bias_name_M1, bias_name_M5, bias_name_M15, bias_name_M30]
    bias_true_all=[bias_true_M1, bias_true_M5, bias_true_M15, bias_true_M30]
    '''
    bias_test_all=[]
    bias_name_all=[]
    bias_true_all=[]
    for m in mass:
        bias_test_all.append(bias_test[m])
        bias_name_all.append(bias_name[m])
        bias_true_all.append(bias_true[m])

    for i in range(len(bias_test_all)):
        bias_test=bias_test_all[i]
        bias_name=bias_name_all[i]
        bias_true=bias_true_all[i]

        x=np.array(bias_test)
        y=0.5*np.ones(len(bias_test))

        plt.xlabel('$\\frac{\mu-\\tilde{\mu}}{\sigma_{\mu}}$')
        plt.ylabel('$\\tilde{\mu}$')
        plt.xlim(xmax=0.5,xmin=-0.5)
        plt.ylim(ymax=1.0,ymin=0)

        colors = ['red','orange','blue','purple','darkturquoise','yellow','green','pink']
        area = np.pi * 4**2

        plt.plot([0.14, 0.14], [0., 1], c='black', linestyle='--')
        plt.plot([-0.14, -0.14], [0., 1], c='black', linestyle='--')
        plt.yticks(np.array([0.5]),np.array(bias_true))
        for j in range(len(bias_test)):
            plt.scatter(x[j], y[j], s=area, c=colors[j], alpha=0.4, label=bias_name[j])

        #plt.title('Toy Function: 3rd Order Power Low (' +mass[i] + ')')
        plt.grid()
        plt.legend()
        plt.savefig('./BiasPlot_UL/'+mass[i] + '.png')
        plt.close('all')
        #plt.show()

main()
