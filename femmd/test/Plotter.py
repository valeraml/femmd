import os
import shutil
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import time

#class FileType:
#    self.XYZ = 0
#    self.VTK = 1

def long_edges(x, y, triangles, radio=4):
    out = []
    for points in triangles:
        #print points
        a,b,c = points
        d0 = np.sqrt( (x[a] - x[b]) **2 + (y[a] - y[b])**2 )
        d1 = np.sqrt( (x[b] - x[c]) **2 + (y[b] - y[c])**2 )
        d2 = np.sqrt( (x[c] - x[a]) **2 + (y[c] - y[a])**2 )
        max_edge = max([d0, d1, d2])
        #print points, max_edge
        if max_edge > radio:
            out.append(True)
        else:
            out.append(False)
    return out

class Plotter:
    def __init__(self, r, v, tris, xmax, ymax, type, sigma=1,**kwargs):
        self.r = r
        self.v = v
        self.tris = tris
        self.dy = 0.1*ymax
        self.dx = 0.1*xmax
        self.xmin = -self.dx
        self.ymin = -self.dy
        self.xmax = xmax+self.dx
        self.ymax = ymax+self.dy
        self.type = type
        self.size = 15
        self.frame_count = 0
        self.savepngs = 0
        self.radius=sigma/2
        if os.path.exists('pngs'):
            print(os.listdir('pngs'))
            #os.system('rm -rf pngs')
            #os.system('del pngs /s')
            shutil.rmtree('pngs')
        os.mkdir("pngs")
        #print(os.listdir('pngs'))

        plt.close()
        #self.fig, self.ax1 = plt.subplots()
        self.fig = plt.figure()
        self.ax1 = self.fig.gca()        
        self.fig.set_size_inches(self.size, self.size)
        self.ax1.set_xlim(self.xmin,self.xmax)
        self.ax1.set_ylim(self.ymin,self.ymax)
        self.ax1.set_aspect('equal')
        vv = (self.v*self.v).sum(1)
        colors = vv[tris].mean(axis=1)
        print(colors.shape, tris.shape)
           
    def update(self,steps,walls=None,title=None):
        if self.type == "pts":
            #self.line1.set_xdata(self.r[:,0])
            #self.line1.set_ydata(self.r[:,1])
            self.ax1.clear()
            self.ax1.set_aspect('equal')
            self.fig.set_size_inches([self.size, self.size],forward=True) 
            self.ax1.set_xlim(self.xmin,self.xmax)
            self.ax1.set_ylim(self.ymin,self.ymax)
            #self.line1, = self.ax1.plot(self.r[:,0],self.r[:,1],'ro',ms=self.screen_dia)
            #self.ax1.scatter(self.r[:,0],self.r[:,1],s=self.particle_size)
            for x, y in zip(self.r[:,0],self.r[:,1]): self.ax1.add_artist(Circle(xy=(x,y),radius=self.radius))
            self.ax1.set_title('steps=%d' % steps)
            plt.pause(0.001)
            pngname = "pngs//frame_%06d.png"% (self.frame_count)
            self.frame_count += 1
            plt.savefig(pngname, dpi=self.fig.dpi);
        elif self.type == "tri":
            vv = (self.v*self.v).sum(1)
            colors = vv[self.tris].mean(axis=1)
            max = colors.max()
            colors *= 1/max;
            self.ax1.clear()
            self.fig.set_size_inches([self.size, self.size]) 
            self.ax1.set_xlim(self.xmin,self.xmax)
            self.ax1.set_ylim(self.ymin,self.ymax)
            #self.ax1.triplot(self.r[:,0],self.r[:,1],self.tris)
            self.ax1.tripcolor(self.r[:,0],self.r[:,1],self.tris,facecolors=colors)
            self.ax1.set_title('steps=%d' % steps)
            plt.pause(0.001)
            pngname = "pngs//frame_%06d.png"% (self.frame_count)
            self.frame_count += 1
            plt.savefig(pngname, dpi=self.fig.dpi);
        elif self.type == "tripts":
            vv = (self.v*self.v).sum(1)
            vv *= 1.0/vv.max()
            colors = vv[self.tris].mean(axis=1)
            max = colors.max()
            colors *= 1/max;
            self.ax1.clear() 
            self.ax1.set_aspect('equal')
            self.fig.set_size_inches([self.size, self.size],forward=True) 
            self.ax1.set_xlim(self.xmin,self.xmax)
            self.ax1.set_ylim(self.ymin,self.ymax)
            self.ax1.axhline(0,self.xmin+self.dx,self.xmax-self.dx, lw = 10)
            self.ax1.axhline(self.ymax - self.dy,self.xmin + self.dx ,self.xmax - self.dx, lw = 10)
            self.ax1.triplot(self.r[:,0],self.r[:,1],self.tris)
            #im=self.ax1.tripcolor(self.r[:,0],self.r[:,1],self.tris,facecolors=colors)
            for x, y, c in zip(self.r[:,0],self.r[:,1], vv): self.ax1.add_artist(Circle(xy=(x,y),fill=False,lw=2,radius=self.radius, ec=(c*.9+.1,0,0)))
            #self.fig.colorbar(im,self.ax1)
            self.ax1.set_title('steps=%d' % steps)
            plt.pause(0.001)
            pngname = "pngs//frame_%06d.png"% (self.frame_count)
            self.frame_count += 1
            plt.savefig(pngname, dpi=self.fig.dpi);
        elif self.type == "tripts1":
            vv = (self.v*self.v).sum(1)
            vv *= 1.0/vv.max()
            colors = vv[self.tris].mean(axis=1)
            max = colors.max()
            colors *= 1/max;
            if walls is None:
                left = self.xmin + self.dx
                right = self.xmax - self.dx
                bottom = self.ymin + self.dy
                top  = self.ymax - self.dy
            else:
                left = walls[0]
                right = walls[1]
                bottom = walls[2]
                top = walls[3]

            self.ax1.clear() 
            self.ax1.set_aspect('equal')
            self.fig.set_size_inches([self.size, self.size],forward=True) 
            self.ax1.set_xlim(self.xmin,self.xmax)
            self.ax1.set_ylim(self.ymin,self.ymax)
            self.ax1.axvline(left,   self.xmin + self.dx, self.xmax - self.dx, lw = 10)
            self.ax1.axvline(right,  self.xmin + self.dx ,self.xmax - self.dx, lw = 10)
            self.ax1.axhline(bottom, self.ymin + self.dy, self.ymax - self.dy, lw = 10)
            self.ax1.axhline(top,    self.ymin + self.dy, self.ymax - self.dy, lw = 10)
            x = self.r[:,0]
            y = self.r[:,1]
            t = self.tris
            self.ax1.triplot(self.r[:,0],self.r[:,1],self.tris,mask=long_edges(x,y,t))
            #self.ax1.tripcolor(self.r[:,0],self.r[:,1],self.tris,facecolors=colors,mask=long_edges(x,y,t))
            #for x, y, c in zip(self.r[:,0],self.r[:,1], vv): self.ax1.add_artist(Circle(xy=(x,y),fill=False,lw=2,radius=self.radius, ec=(c*.9+.1,0,0)))
            if title is None:
                self.ax1.set_title('steps=%d' % steps)
            else:
               self.ax1.set_title(title) 
            self.fig.canvas.draw()
            #self.fig.canvas.flush_events()
            plt.pause(0.001)
            pngname = "pngs//frame_%06d.png"% (self.frame_count)
            self.frame_count += 1
            plt.savefig(pngname, dpi=self.fig.dpi);
        else:
            print("unknown option %s\n" % self.type)
            exit()
        #time.sleep(0.01)
