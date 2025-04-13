from ROOT import TH1F, TH1, TH2F, TH2
import ROOT
from numpy import sqrt, log, ceil, log2
from tqdm import tqdm
import numpy as np
from pdb import set_trace

def calc_sig(sig, bkg,s_err,b_err):

    ntot = sig + bkg

    #   if(sig <= 0): return 0, 0
    #   if(bkg <= 0): return 0, 0

    #Counting experiment
    significance = sqrt(2*((ntot * log(ntot/bkg)) - sig)) # why a 2 is here
    #   significance = sig / sqrt(bkg)

    #error on significance
    uncert = sqrt((log(ntot/bkg)*s_err)**2 + ((log(ntot/bkg) - (sig/bkg))*b_err)**2)
    #   uncert = sqrt((s_err**2/bkg) + (significance/2/bkg*b_err)**2)

    return significance, uncert

def sum_z(zs):
    sumu=0
    sumz=0
    for i in range(len(zs)):
        sumz+=zs[i][0]**2
        sumu+=(zs[i][0]*zs[i][1])**2
    sumz=sqrt(sumz)
    sumu=sqrt(sumu)/sumz if sumz != 0 else 0
    return sumz, sumu

def copy_hist(hname, h0, bl, br):

    nbin = br - bl + 1
    l_edge = h0.GetBinLowEdge(bl)
    r_edge = l_edge + h0.GetBinWidth(bl) * nbin
    h1 = TH1F(hname, hname, nbin, l_edge, r_edge)

    for i, j in enumerate(range(bl, br+1)):

        y = h0.GetBinContent(j)
        err = h0.GetBinError(j)

        h1.SetBinContent(i+1, y)
        h1.SetBinError(i+1, err)
        
    return h1

class pyTH1(object):
    """Customized TH1 class"""

    def __init__(self, hist):

        assert isinstance(hist, TH1), 'hist==ROOT.TH1F'

        self.BinContent = np.array([])
        self.BinError = np.array([])
        self.BinLowEdge = np.array([])
        self.BinWidth = np.array([])

        for i in range(1, hist.GetSize()-1):

            self.BinContent = np.append(self.BinContent, hist.GetBinContent(i))
            self.BinError = np.append(self.BinError, hist.GetBinError(i))
            self.BinLowEdge = np.append(self.BinLowEdge, hist.GetBinLowEdge(i))
            self.BinWidth = np.append(self.BinWidth, hist.GetBinWidth(i))

    def GetSize(self):

        return len(self.BinContent)

    def GetLowEdge(self):

        return self.BinLowEdge[0]

    def GetHighEdge(self):

        return self.BinLowEdge[-1] + self.BinWidth[-1]

    def GetBinContent(self, i=0):

        if i==0: return self.BinContent
        return self.BinContent[i-1]

    def GetBinError(self, i=0):

        if i==0: return self.BinError
        return self.BinError[i-1]

    def GetBinLowEdge(self, i=0):

        if i==0: return self.BinLowEdge
        return self.BinLowEdge[i-1]

    def GetBinWidth(self, i=0):

        if i==0: return self.BinWidth
        return self.BinWidth[i-1]

    def Integral(self, i=0, j=0):

        if i==0 and j==0:
            return self.BinContent.sum()

        return self.BinContent[i-1:j].sum()

    def IntegralAndError(self, i=0, j=0):

        if i==0 and j==0:
            return self.BinContent.sum(), sqrt((self.BinError**2).sum())

        return self.BinContent[i-1:j].sum(), sqrt((self.BinError[i-1:j]**2).sum())

    def to_TH1F(self, name, bl=-1, br=-1):

        if bl == -1 and br == -1:
            bl = 1
            br = self.GetSize()

        hist = TH1F(name, name, br - bl + 1, self.GetBinLowEdge(bl), self.GetBinLowEdge(br) + self.GetBinWidth(br))

        for i, j in enumerate(range(bl, br+1)):

            y = self.GetBinContent(j)
            err = self.GetBinError(j)

            hist.SetBinContent(i+1, y)
            hist.SetBinError(i+1, err)

        return hist
    
    def Smooth(self, iterations):
        # Create a temporary TH1F histogram to perform smoothing
        temp_hist = ROOT.TH1F("temp_hist", "Temporary Histogram", len(self.BinContent), self.BinLowEdge[0], self.BinLowEdge[-1] + self.BinWidth[-1])
        
        # Fill the temporary histogram with the bin contents
        for i in range(len(self.BinContent)):
            temp_hist.SetBinContent(i+1, self.BinContent[i])

        # Smooth the temporary histogram
        temp_hist.Smooth(iterations)

        # Update the bin contents after smoothing
        for i in range(len(self.BinContent)):
            self.BinContent[i] = temp_hist.GetBinContent(i+1)
            self.BinError[i] = temp_hist.GetBinError(i+1)

class pyTH2(object):
    """Customized TH2 class"""

    def __init__(self, hist):
        assert isinstance(hist, TH2F), 'hist should be an instance of ROOT.TH2F'

        self.BinContent = np.array([])
        self.BinError = np.array([])
        self.BinXLowEdge = np.array([])
        self.BinXWidth = np.array([])
        self.BinYLowEdge = np.array([])
        self.BinYWidth = np.array([])

        for i in range(1, hist.GetNbinsX()+1):
            for j in range(1, hist.GetNbinsY()+1):
                self.BinContent = np.append(self.BinContent, hist.GetBinContent(i, j))
                self.BinError = np.append(self.BinError, hist.GetBinError(i, j))
                if j == 1:
                    self.BinXLowEdge = np.append(self.BinXLowEdge, hist.GetXaxis().GetBinLowEdge(i))
                    self.BinXWidth = np.append(self.BinXWidth, hist.GetXaxis().GetBinWidth(i))
                if i == 1:
                    self.BinYLowEdge = np.append(self.BinYLowEdge, hist.GetYaxis().GetBinLowEdge(j))
                    self.BinYWidth = np.append(self.BinYWidth, hist.GetYaxis().GetBinWidth(j))

    def GetSize(self):
        return len(self.BinContent)

    def GetXSize(self):
        return len(self.BinXLowEdge)

    def GetYSize(self):
        return len(self.BinYLowEdge)

    def GetXLowEdge(self):
        return self.BinXLowEdge[0]

    def GetXHighEdge(self):
        return self.BinXLowEdge[-1] + self.BinXWidth[-1]

    def GetYLowEdge(self):
        return self.BinYLowEdge[0]

    def GetYHighEdge(self):
        return self.BinYLowEdge[-1] + self.BinYWidth[-1]

    def GetBinContent(self, i=0, j=0):
        if i == 0 and j == 0:
            return self.BinContent
        index = (i-1) * len(self.BinYLowEdge) + (j-1)
        return self.BinContent[index]

    def GetBinError(self, i=0, j=0):
        if i == 0 and j == 0:
            return self.BinError
        index = (i-1) * len(self.BinYLowEdge) + (j-1)
        return self.BinError[index]

    def GetXBinLowEdge(self, i=0):
        if i == 0:
            return self.BinXLowEdge
        return self.BinXLowEdge[i-1]

    def GetXBinWidth(self, i=0):
        if i == 0:
            return self.BinXWidth
        return self.BinXWidth[i-1]

    def GetYBinLowEdge(self, i=0):
        if i == 0:
            return self.BinYLowEdge
        return self.BinYLowEdge[i-1]

    def GetYBinWidth(self, i=0):
        if i == 0:
            return self.BinYWidth
        return self.BinYWidth[i-1]

    def Integral(self, ix1=0, ix2=0, iy1=0, iy2=0):
        if ix1 == 0 and ix2 == 0 and iy1 == 0 and iy2 == 0:
            return self.BinContent.sum()

        ix1 = max(ix1 - 1, 0)
        ix2 = min(ix2, len(self.BinXLowEdge))
        iy1 = max(iy1 - 1, 0)
        iy2 = min(iy2, len(self.BinYLowEdge))

        total = 0
        for i in range(ix1, ix2):
            for j in range(iy1, iy2):
                total += self.BinContent[i * len(self.BinYLowEdge) + j]
        return total
    
    def IntegralX(self, ix1=0, ix2=0):
        if ix1 == 0 and ix2 == 0:
            ix1, ix2 = 1, len(self.BinXLowEdge)
        ix1 = max(ix1 - 1, 0)
        ix2 = min(ix2, len(self.BinXLowEdge))
        total = np.zeros(len(self.BinYLowEdge))
        for j in range(len(self.BinYLowEdge)):
            index = ix1 + j * len(self.BinXLowEdge)
            total[j] += self.BinContent[index:index+ix2-ix1].sum()
        return total

    def IntegralAndError(self, ix1=0, ix2=0, iy1=0, iy2=0):
        if ix1 == 0 and ix2 == 0 and iy1 == 0 and iy2 == 0:
            return self.BinContent.sum(), sqrt((self.BinError**2).sum())

        ix1 = max(ix1 - 1, 0)
        ix2 = min(ix2, len(self.BinXLowEdge))
        iy1 = max(iy1 - 1, 0)
        iy2 = min(iy2, len(self.BinYLowEdge))

        total = 0
        error2 = 0
        for i in range(ix1, ix2):
            for j in range(iy1, iy2):
                index = i * len(self.BinYLowEdge) + j
                total += self.BinContent[index]
                error2 += self.BinError[index]**2
        return total, sqrt(error2)

    def IntegralXAndError(self, ix1=0, ix2=0):
        if ix1 == 0 and ix2 == 0:
            ix1, ix2 = 1, len(self.BinXLowEdge)
        ix1 = max(ix1 - 1, 0)
        ix2 = min(ix2, len(self.BinXLowEdge))
        total = np.zeros(len(self.BinYLowEdge))
        error2 = np.zeros(len(self.BinYLowEdge))
        for j in range(len(self.BinYLowEdge)):
            index = ix1 + j * len(self.BinXLowEdge)
            total[j] = self.BinContent[index:index+ix2-ix1].sum()
            error2[j] = (self.BinError[index:index+ix2-ix1]**2).sum()
        
        return total, sqrt(error2)

    def to_TH2F(self, name, xl=-1, xr=-1, yl=-1, yr=-1):
        if xl == -1 and xr == -1:
            xl = 1
            xr = len(self.BinXLowEdge)
        if yl == -1 and yr == -1:
            yl = 1
            yr = len(self.BinYLowEdge)

        hist = TH2F(name, name, xr - xl + 1, self.GetXBinLowEdge(xl), self.GetXBinLowEdge(xr) + self.GetXBinWidth(xr), yr - yl + 1, self.GetYBinLowEdge(yl), self.GetYBinLowEdge(yr) + self.GetYBinWidth(yr))

        for ix in range(xl, xr + 1):
            for iy in range(yl, yr + 1):
                index = (ix - 1) * len(self.BinYLowEdge) + (iy - 1)
                y = self.GetBinContent(ix, iy)
                err = self.GetBinError(ix, iy)
                hist.SetBinContent(ix - xl + 1, iy - yl + 1, y)
                hist.SetBinError(ix - xl + 1, iy - yl + 1, err)

        return hist
    
    def Smooth(self, iterations):
        # Create a temporary TH2F histogram to perform smoothing
        temp_hist = ROOT.TH2F("temp_hist", "Temporary Histogram", len(self.BinXLowEdge), self.BinXLowEdge[0], self.BinXLowEdge[-1] + self.BinXWidth[-1], len(self.BinYLowEdge), self.BinYLowEdge[0], self.BinYLowEdge[-1] + self.BinYWidth[-1])
        
        # Fill the temporary histogram with the bin contents
        for ix in range(len(self.BinXLowEdge)):
            for iy in range(len(self.BinYLowEdge)):
                index = ix * len(self.BinYLowEdge) + iy
                temp_hist.SetBinContent(ix+1, iy+1, self.BinContent[index])

        # Smooth the temporary histogram
        temp_hist.Smooth(iterations)

        # Update the bin contents after smoothing
        for ix in range(len(self.BinXLowEdge)):
            for iy in range(len(self.BinYLowEdge)):
                index = ix * len(self.BinYLowEdge) + iy
                self.BinContent[index] = temp_hist.GetBinContent(ix+1, iy+1)
                self.BinError[index] = temp_hist.GetBinError(ix+1, iy+1)

class categorizer(object):

    def __init__(self, h_sig, h_bkg, h_bkg_rw_num=None, h_bkg_rw_den=None):

        assert isinstance(h_sig, TH1), 'h_sig==ROOT.TH1'
        assert isinstance(h_bkg, TH1), 'h_bkg==ROOT.TH1'
        assert h_sig.GetSize() == h_bkg.GetSize(), 'h_sig and h_bkg have different sizes'
        assert h_sig.GetBinLowEdge(1) == h_bkg.GetBinLowEdge(1), 'h_sig and h_bkg have different edges'
        assert h_sig.GetBinWidth(1) == h_bkg.GetBinWidth(1), 'h_sig and h_bkg have different bin widths'
        self.reweight = False
        if h_bkg_rw_num and h_bkg_rw_den:
            self.reweight = True
            assert h_sig.GetSize() == h_bkg_rw_num.GetSize(), 'h_sig and h_bkg have different sizes'
            assert h_sig.GetBinLowEdge(1) == h_bkg_rw_num.GetBinLowEdge(1), 'h_sig and h_bkg have different edges'
            assert h_sig.GetBinWidth(1) == h_bkg_rw_num.GetBinWidth(1), 'h_sig and h_bkg have different bin widths'
            assert h_sig.GetSize() == h_bkg_rw_den.GetSize(), 'h_sig and h_bkg have different sizes'
            assert h_sig.GetBinLowEdge(1) == h_bkg_rw_den.GetBinLowEdge(1), 'h_sig and h_bkg have different edges'
            assert h_sig.GetBinWidth(1) == h_bkg_rw_den.GetBinWidth(1), 'h_sig and h_bkg have different bin widths'

        self.h_sig = pyTH1(h_sig)
        self.h_bkg = pyTH1(h_bkg)
        if self.reweight:
            self.h_bkg_rw_num = pyTH1(h_bkg_rw_num)
            self.h_bkg_rw_den = pyTH1(h_bkg_rw_den)

    def smooth_sim(self, iterations, SorB='B'):

        if SorB == 'S': py_hist = self.h_sig
        elif SorB == 'B': py_hist = self.h_bkg

        hist = py_hist.to_TH1F("hist")

        py_hist.Smooth(iterations)

        if SorB == 'S': self.h_sig = py_hist
        elif SorB == 'B': self.h_bkg = py_hist

        # Create canvas to draw histograms
        canvas = ROOT.TCanvas("canvas", "Smoothed Histogram", 800, 600)

        # Draw original histogram
        hist.SetLineColor(ROOT.kBlue)
        hist.Draw()

        # Set y-axis to logarithmic scale
        canvas.SetLogy(True)

        # Draw smoothed histogram
        hist_smoothed = ROOT.TH1F("hist_smoothed", "Smoothed Histogram", len(py_hist.BinContent), py_hist.BinLowEdge[0], py_hist.BinLowEdge[-1] + py_hist.BinWidth[-1])
        for i in range(len(py_hist.BinContent)):
            hist_smoothed.SetBinContent(i+1, py_hist.BinContent[i])
        hist_smoothed.SetLineColor(ROOT.kRed)
        hist_smoothed.Draw("same")

        # Create legend
        legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
        legend.AddEntry(hist, "Original", "l")
        legend.AddEntry(hist_smoothed, "Smoothed", "l")
        legend.Draw()

        # Update canvas
        canvas.Modified()
        canvas.Update()

        # Save canvas as PNG file
        canvas.SaveAs("smoothed_histogram.png")

    def smooth(self, bl, br, SorB='B', function='Epoly2'):

        if SorB == 'S': hist = self.h_sig
        elif SorB == 'B': hist = self.h_bkg

        nbin = hist.GetSize()
        l_edge = hist.GetLowEdge()
        r_edge = hist.GetHighEdge()

        h_merge_list = ROOT.TList()

        if bl == 0:
            print("!!! Left edge can not be 0 !!!")
            return

        if bl > 1:
            h_left = hist.to_TH1F('h_left', 1, bl-1)
            h_merge_list.Add(h_left)

        hist_to_smooth = hist.to_TH1F('hist_to_smooth', bl, br)
        smoothed_hist = fit_BDT('smoothed_hist', hist_to_smooth, function)
        h_merge_list.Add(smoothed_hist)

        if br != nbin:
            h_right = hist.to_TH1F('h_right', br+1, nbin)
            h_merge_list.Add(h_right)

        htemp = TH1F('htemp', 'htemp', nbin, l_edge, r_edge)
        htemp.Merge(h_merge_list)

        if SorB == 'S':
            self.h_sig = pyTH1(htemp)
        elif SorB == 'B':
            self.h_bkg = pyTH1(htemp)

        if bl != 1: h_left.Delete()
        hist_to_smooth.Delete()
        smoothed_hist.Delete()
        if br != nbin: h_right.Delete()
        htemp.Delete()

    def fit(self, bl, br, nbin, minN=5, floatB=False, earlystop=-1, pbar=False):

        if nbin == 1:

            if floatB: return [], 0, 0

            nsig, dsig = self.h_sig.IntegralAndError(bl, br)
            nbkg, dbkg = self.h_bkg.IntegralAndError(bl, br)
            # temp_nbkg = nbkg
            if self.reweight and nbkg != 0:
                if self.h_bkg_rw_den.Integral(bl, br) == 0: print("what!!!", bl, br)
                if self.h_bkg_rw_den.Integral(bl, br) < 5*minN or self.h_bkg_rw_num.Integral(bl, br) < 5*minN: return -1, -1, -1
                nbkg, dbkg = nbkg*self.h_bkg_rw_num.Integral(bl, br)/self.h_bkg_rw_den.Integral(bl, br), dbkg*self.h_bkg_rw_num.Integral(bl, br)/self.h_bkg_rw_den.Integral(bl, br)
                
            if nbkg < minN: return -1, -1, -1

            # print(bl, br, nbkg, temp_nbkg, self.h_bkg_rw_num.Integral(bl, br), self.h_bkg_rw_den.Integral(bl, br), sep=", ")

            z, u = calc_sig(nsig, nbkg, dsig, dbkg)

            return [bl], z, u

        elif nbin > 1:

            L = int(ceil(log2(nbin)))
            N2 = 2**(L - 1)
            N1 = nbin - N2

            bmax, zmax, umax, stop = -1, -1, -1, 0
 
            for b in (range(bl, br+1) if not pbar else tqdm(range(bl, br+1), ncols=70)):

                b1, z1, u1 = self.fit(bl, b-1, N1, minN=minN, floatB=floatB, earlystop=earlystop)
                if b1 == -1: continue
                b2, z2, u2 = self.fit(b, br, N2, minN=minN, earlystop=earlystop)
                if b2 == -1: break

                z = sqrt(z1**2 + z2**2)
                u = sqrt((u1*z1)**2 + (u2*z2)**2)/z if z != 0 else 0

                if z > zmax:
                    stop = 0
                    zmax = z
                    bmax = sorted(list(set(b1 + [b] + b2)))
                    umax = u
                else:
                    stop += 1
                if stop == earlystop: break

            return bmax, zmax, umax

class categorizer_shape(object):

    def __init__(self, h_sig, h_bkg, h_bkg_rw_num=None, h_bkg_rw_den=None):
        assert isinstance(h_sig, TH2), 'h_sig should be an instance of TH2'
        assert isinstance(h_bkg, TH2), 'h_bkg should be an instance of TH2'
        assert h_sig.GetXaxis().GetNbins() == h_bkg.GetXaxis().GetNbins(), 'h_sig and h_bkg have different sizes'
        assert h_sig.GetXaxis().GetBinLowEdge(1) == h_bkg.GetXaxis().GetBinLowEdge(1), 'h_sig and h_bkg have different edges'
        assert h_sig.GetXaxis().GetBinWidth(1) == h_bkg.GetXaxis().GetBinWidth(1), 'h_sig and h_bkg have different bin widths'
        self.reweight = False
        if h_bkg_rw_num and h_bkg_rw_den:
            self.reweight = True
            assert h_sig.GetXaxis().GetNbins() == h_bkg_rw_num.GetXaxis().GetNbins(), 'h_sig and h_bkg have different sizes'
            assert h_sig.GetXaxis().GetBinLowEdge(1) == h_bkg_rw_num.GetBinLowEdge(1), 'h_sig and h_bkg have different edges'
            assert h_sig.GetXaxis().GetBinWidth(1) == h_bkg_rw_num.GetBinWidth(1), 'h_sig and h_bkg have different bin widths'
            assert h_sig.GetXaxis().GetNbins() == h_bkg_rw_den.GetXaxis().GetNbins(), 'h_sig and h_bkg have different sizes'
            assert h_sig.GetXaxis().GetBinLowEdge(1) == h_bkg_rw_den.GetBinLowEdge(1), 'h_sig and h_bkg have different edges'
            assert h_sig.GetXaxis().GetBinWidth(1) == h_bkg_rw_den.GetBinWidth(1), 'h_sig and h_bkg have different bin widths'

        self.h_sig = pyTH2(h_sig)
        self.h_bkg = pyTH2(h_bkg)
        if self.reweight:
            self.h_bkg_rw_num = pyTH1(h_bkg_rw_num)
            self.h_bkg_rw_den = pyTH1(h_bkg_rw_den)

    def smooth_sim(self, iterations, SorB='B'):
        if SorB == 'S': py_hist = self.h_sig
        elif SorB == 'B': py_hist = self.h_bkg

        hist = py_hist.to_TH2F("hist")

        py_hist.Smooth(iterations)

        if SorB == 'S': self.h_sig = py_hist
        elif SorB == 'B': self.h_bkg = py_hist

        # Create canvas to draw histograms
        canvas = ROOT.TCanvas("canvas", "Smoothed Histogram", 800, 600)

        # Draw original histogram
        hist.SetLineColor(ROOT.kBlue)
        hist.Draw("COLZ")

        # Draw smoothed histogram
        hist_smoothed = ROOT.TH2F("hist_smoothed", "Smoothed Histogram", len(py_hist.BinXLowEdge), py_hist.BinXLowEdge[0], py_hist.BinXLowEdge[-1] + py_hist.BinXWidth[-1], len(py_hist.BinYLowEdge), py_hist.BinYLowEdge[0], py_hist.BinYLowEdge[-1] + py_hist.BinYWidth[-1])
        for ix in range(len(py_hist.BinXLowEdge)):
            for iy in range(len(py_hist.BinYLowEdge)):
                hist_smoothed.SetBinContent(ix + 1, iy + 1, py_hist.BinContent[ix * len(py_hist.BinYLowEdge) + iy])
        hist_smoothed.SetLineColor(ROOT.kRed)
        hist_smoothed.Draw("same COLZ")

        # Create legend
        legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
        legend.AddEntry(hist, "Original", "l")
        legend.AddEntry(hist_smoothed, "Smoothed", "l")
        legend.Draw()

        # Update canvas
        canvas.Modified()
        canvas.Update()

        # Save canvas as PNG file
        canvas.SaveAs("smoothed_histogram.png")

    def smooth(self, bl, br, SorB='B', function='Epoly2'):
        if SorB == 'S': hist = self.h_sig
        elif SorB == 'B': hist = self.h_bkg

        nbin = hist.GetSize()
        l_edge = hist.GetXLowEdge()
        r_edge = hist.GetXHighEdge()

        h_merge_list = ROOT.TList()

        if bl == 0:
            print("!!! Left edge can not be 0 !!!")
            return

        if bl > 1:
            h_left = hist.to_TH2F('h_left', 1, bl-1)
            h_merge_list.Add(h_left)

        hist_to_smooth = hist.to_TH2F('hist_to_smooth', bl, br)
        smoothed_hist = fit_BDT('smoothed_hist', hist_to_smooth, function)
        h_merge_list.Add(smoothed_hist)

        if br != nbin:
            h_right = hist.to_TH2F('h_right', br+1, nbin)
            h_merge_list.Add(h_right)

        htemp = TH2F('htemp', 'htemp', nbin, l_edge, r_edge, len(hist.BinYLowEdge), hist.BinYLowEdge[0], hist.BinYLowEdge[-1] + hist.BinYWidth[-1])
        htemp.Merge(h_merge_list)

        if SorB == 'S':
            self.h_sig = pyTH2(htemp)
        elif SorB == 'B':
            self.h_bkg = pyTH2(htemp)

        if bl != 1: h_left.Delete()
        hist_to_smooth.Delete()
        smoothed_hist.Delete()
        if br != nbin: h_right.Delete()
        htemp.Delete()

    def fit(self, bl, br, nbin, minN=5, floatB=False, earlystop=-1, pbar=False):
        if nbin == 1:
            if floatB: return [], 0, 0

            nsig, dsig = self.h_sig.IntegralXAndError(bl, br)
            nbkg, dbkg = self.h_bkg.IntegralXAndError(bl, br)
            if self.reweight and nbkg.sum() != 0:
                if self.h_bkg_rw_den.Integral(bl, br) == 0: print("what!!!", bl, br)
                nbkg, dbkg = nbkg * self.h_bkg_rw_num.Integral(bl, br) / self.h_bkg_rw_den.Integral(bl, br), dbkg * self.h_bkg_rw_num.Integral(bl, br) / self.h_bkg_rw_den.Integral(bl, br)

            if nbkg.sum() < minN: return -1, -1, -1

            zs, us = calc_sig(nsig, nbkg, dsig, dbkg)

            # z, u = np.sqrt(zs**2), sqrt((us**2 * zs**2).sum()) / np.sqrt(zs**2).sum()
            z, u = zs.mean()*np.sqrt(len(zs)), sqrt((us**2).sum()/len(zs))

            return [bl], z, u

        elif nbin > 1:
            L = int(ceil(log2(nbin)))
            N2 = 2 ** (L - 1)
            N1 = nbin - N2

            bmax, zmax, umax, stop = -1, -1, -1, 0

            for b in (range(bl, br + 1) if not pbar else tqdm(range(bl, br + 1), ncols=70)):
                b1, z1, u1 = self.fit(bl, b-1, N1, minN=minN, floatB=floatB, earlystop=earlystop)
                if b1 == -1: continue
                b2, z2, u2 = self.fit(b, br, N2, minN=minN, earlystop=earlystop)
                if b2 == -1: break

                z = sqrt(z1 ** 2 + z2 ** 2)
                u = sqrt((u1 * z1)**2 + (u2 * z2)**2)

                if z > zmax:
                    stop = 0
                    zmax = z
                    bmax = sorted(list(set(b1 + [b] + b2)))
                    umax = u
                else:
                    stop += 1
                if stop == earlystop: break

            return bmax, zmax, umax

def fit_BDT(hname, hist, function='Epoly2', printMessage=False):

    if not printMessage: ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)  #WARNING

    n_events = hist.Integral()

    nbin = hist.GetSize() - 2
    l_edge = hist.GetBinLowEdge(1)
    r_edge = l_edge + hist.GetBinWidth(1) * nbin

    bdt = ROOT.RooRealVar("bdt", "bdt", l_edge, r_edge)

    a = ROOT.RooRealVar("a", "a", -100, 100)
    b = ROOT.RooRealVar("b", "b", -100, 100)
    c = ROOT.RooRealVar("c", "c", -100, 100)
    d = ROOT.RooRealVar("d", "d", -100, 100)
    e = ROOT.RooRealVar("e", "e", -100, 100)
    f = ROOT.RooRealVar("f", "f", -100, 100)

    pdf = {}
    pdf['power'] = ROOT.RooGenericPdf("pdf","PDF for BDT distribution", "@0**@1", ROOT.RooArgList(bdt, a)) #power

    # EPoly Family
    pdf['Exp'] = ROOT.RooGenericPdf("pdf","PDF for BDT distribution", "exp(@1*@0)", ROOT.RooArgList(bdt, a))
    pdf['Epoly2'] = ROOT.RooGenericPdf("pdf","PDF for BDT distribution", "exp(@1*@0+@2*@0**2)", ROOT.RooArgList(bdt, a, b))
    pdf['Epoly3'] = ROOT.RooGenericPdf("pdf","PDF for BDT distribution", "exp(@1*@0+@2*@0**2+@3*@0**3)", ROOT.RooArgList(bdt, a, b, c))
    pdf['Epoly4'] = ROOT.RooGenericPdf("pdf","PDF for BDT distribution", "exp(@1*@0+@2*@0**2+@3*@0**3+@4*@0**4)", ROOT.RooArgList(bdt, a, b, c, d))
    pdf['Epoly5'] = ROOT.RooGenericPdf("pdf","PDF for BDT distribution", "exp(@1*@0+@2*@0**2+@3*@0**3+@4*@0**4+@5*@0**5)", ROOT.RooArgList(bdt, a, b, c, d, e))
    pdf['Epoly6'] = ROOT.RooGenericPdf("pdf","PDF for BDT distribution", "exp(@1*@0+@2*@0**2+@3*@0**3+@4*@0**4+@5*@0**5+@6*@0**6)", ROOT.RooArgList(bdt, a, b, c, d, e, f))

    data = ROOT.RooDataHist("dh", "dh", ROOT.RooArgList(bdt), hist)
    #data.Print("all")

    pdf[function].fitTo(data, ROOT.RooFit.Verbose(False), ROOT.RooFit.PrintLevel(-1))

    if printMessage:
        a.Print()
        b.Print()
        c.Print()
        d.Print()
        e.Print()
        f.Print()
    
        frame = bdt.frame()
        data.plotOn(frame)
        pdf[function].plotOn(frame)
    
        frame.Draw()
        #frame.Print('test.pdf')
    
        dof = {'power': 2, 'Exp': 2, 'Epoly2': 3, 'Epoly3': 4, 'Epoly4': 5, 'Epoly5': 6, 'Epoly6': 7}
        reduced_chi_square = frame.chiSquare(dof[function])
        probability = TMath.Prob(frame.chiSquare(dof[function]) * (nbin - dof[function]), nbin - dof[function])
        print('chi square:', reduced_chi_square)
        print('probability: ', probability)

    #raw_input('Press enter to continue')

    # fill the fitted pdf into a histogram
    hfit = TH1F(hname, hname, nbin, l_edge, r_edge)
    pdf[function].fillHistogram(hfit, ROOT.RooArgList(bdt), n_events)

    return hfit
