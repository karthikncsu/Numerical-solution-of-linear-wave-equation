"""Solution of one dimensional linear wave equation: ut+a*ux=0 using different
integration schemes. This program aids as a starter for mechanical engineers
who are passionate about numerical programming.

The classes and objects of python are used in this program and the OOPs concepts
help in writing a compact program. The solution of the linear wave equation is
obtained using the following schemes
    
    * First-order Backward Difference Explicit Integration (FOBD_explicit)
    * Second-order Backward Difference Explicit Integration (SOBD_explicit)
    * Crank-Nicolson, Implicit Integration (CN_implicit)
    * Lax-Wendroff Explicit Integration (LW_explicit)
    * Lax Explicit Integration (LAX_explicit)
    * First-order Backward Difference, Implicit Integration (FOBD_implicit)

Separate classes were written for each scheme and the name of the class for  is
as shown in the brackets following each scheme. All the classes have a similar structure
with a constructor and a function with the name evaluate. Common class with the name 
"LinearWaveEqn" is used to perform numerical integration for all schemes.  
"""
import numpy as np
import matplotlib.pyplot as plt
import copy

class FOBD_explicit:
    """
    Class: First-order Backward Difference Explicit Integration

    Args:
        CFL0: CFL number
    
    Methods: evaluate()
    """
    def __init__(self,CFL0):
        self.a= CFL0; self.b=1-CFL0
    def evaluate(self,N0,usol,usoln):
        for i in range(1,N0-2):
            usoln[i,0]=self.a*usol[i-1,0]+self.b*usol[i,0]
        usoln[0,0]=self.b*usol[0,0]
        return usoln

class SOBD_explicit:
    """
    Class: Second-order Backward Difference Explicit Integration

    Args:
        CFL0: CFL number
    
    Methods: evaluate()
    """
    def __init__(self,CFL0):
        self.a=-0.5*CFL0; self.b=2*CFL0; self.c=1-1.5*CFL0;
    def evaluate(self,N0,usol,usoln):
        for i in range(2,N0-2):
            usoln[i,0]=self.a*usol[i-2,0]+self.b*usol[i-1,0]+self.c*usol[i,0];
        usoln[0,0]=(self.a+self.c)*usol[0,0];
        usoln[1,0]=self.b*usol[0,0]+self.c*usol[1,0];
        return usoln

class CN_implicit:
    """
    Class: Crank-Nicolson, Implicit Integration

    Args:
        CFL0: CFL number
    
    Methods: evaluate()
    """
    def __init__(self,CFL0,N0):
        K=np.zeros((N0-2,N0-2))
        K1=np.zeros((N0-2,N0-2))
        for i in range(1,N0-3):
            K[i,i-1]=-CFL0/4.0;  K[i,i]=1.0;    K[i,i+1]=CFL0/4.0;
            K1[i,i-1]=CFL0/4.0;  K1[i,i]=1.0;   K1[i,i+1]=-CFL0/4.0;
        K[0,0]=1.0; K[0,1]=CFL0/4.0;  K[N0-3,N0-4]=-CFL0/4.0;  K[N0-3,N0-3]=1.0
        K1[0,0]=1.0; K1[0,1]=-CFL0/4.0;  K1[N0-3,N0-4]=CFL0/4.0;  K1[N0-3,N0-3]=1

        self.K=K; self.K1=K1
    def evaluate(self,N0,usol,usoln):
        f=np.matmul(self.K1,usol);
        usoln=thomas_algorithm(self.K,f);
        return usoln

class LW_explicit:
    """
    Class:  Lax-Wendroff Explicit Integration

    Args:
        CFL0: CFL number
    
    Methods: evaluate()
    """
    def __init__(self,CFL0):
        self.a=0.5*CFL0*(1+CFL0); self.b=1-CFL0*CFL0; self.c=0.5*CFL0*(CFL0-1)
    def evaluate(self,N0,usol,usoln):
        for i in range(1,N0-3):
            usoln[i,0]=self.a*usol[i-1,0]+self.b*usol[i,0]+self.c*usol[i+1,0];
        usoln[0,0]=self.b*usol[0,0]+self.c*usol[1,0]
        usoln[N0-3,0]=self.a*usol[N0-4,0]+self.b*usol[N0-3,0];
        return usoln

class LAX_explicit:
    """
    Class: Lax Explicit Integration

    Args:
        CFL0: CFL number
    
    Methods: evaluate()
    """
    def __init__(self,CFL0):
        self.a=0.5*(1+CFL0); self.b=0.5*(1-CFL0);
    def evaluate(self,N0,usol,usoln):
        for i in range(1,N0-3):
            usoln[i,0]=self.a*usol[i-1,0]+self.b*usol[i+1,0];
        usoln[0,0]=self.b*usol[1,0]
        usoln[N0-3,0]=self.b*usol[N0-3,0];
        return usoln

class FOBD_implicit:
    """
    Class: First-order Backward Difference, Implicit Integration

    Args:
        CFL0: CFL number
    
    Methods: evaluate()
    """
    def __init__(self,CFL0):
        self.a=-CFL0; self.b=1+CFL0
    def evaluate(self,N0,usol,usoln):
        usoln[0,0]=usol[0,0]/self.b;
        for i in range(1,N0-3):
            usoln[i,0]=(usol[i,0]-self.a*usoln[i-1,0])/self.b;
        return usoln

def thomas_algorithm(K,f):
    """
    Function: Thomas Algorithm to Solve Tridiagonal System of Matrices

    Args:
        K,f: System matrices
    """
    K=copy.deepcopy(K)
    f=copy.deepcopy(f)
    (N1,N2)=K.shape
    usol=np.zeros((N1,1))
    # Eliminating right Diagnol Elements
    for i in range(1,N1):
        K[i,i]=K[i,i]-K[i-1,i]/K[i-1,i-1]*K[i,i-1]
        f[i,0]=f[i,0]-f[i-1,0]/K[i-1,i-1]*K[i,i-1]
        K[i,i-1]=0.0;

    usol[N1-1,0]=f[N1-1,0]/K[N1-1,N1-1]
    for i in np.arange(N1-2,-1,-1):
        usol[i,0]=(f[i,0]-usol[i+1,0]*K[i,i+1])/K[i,i]
    return usol

class LinearWaveEqn:
    """
    Class: Solution of wave equation using different schemes

    Methods: assemblematrix(), solution()
    """
    def assemblematrix(self,scheme,CFL0,N0):
        """
        Function: Defines objects for different schemes based on the scheme selected

        Args:
            method: Integration scheme
        """
        if scheme=="FOBD_explicit":
            self.solnscheme=FOBD_explicit(CFL0)
            figtitle='First-order Backward Difference Explicit Integration '
            figname='FOBD_explicit_integration_'
        elif scheme=="SOBD_explicit":
            self.solnscheme=SOBD_explicit(CFL0)
            figtitle='Second-order Backward Difference Explicit Integration '
            figname='SOBD_explicit_integration_'
        elif scheme=="CN_implicit":
            self.solnscheme=CN_implicit(CFL0,N0)
            figtitle='Crank-Nicolson, Implicit Integration '
            figname='Crank_Nicolson_implicit_integration_'
        elif scheme=="LW_explicit":
            self.solnscheme=LW_explicit(CFL0)
            figtitle='Lax-Wendroff Explicit Integration Scheme '
            figname='Laxwendroff_imlicit_integration_'
        elif scheme=="LAX_explicit":
            self.solnscheme=LAX_explicit(CFL0)
            figtitle='Lax Explicit Integration Scheme '
            figname='lax_scheme_integration_'
        elif scheme=="FOBD_implicit":
            self.solnscheme=FOBD_implicit(CFL0)
            figtitle='First-order Backward Difference, Implicit Integration Scheme '
            figname='FOBD_implicit_integration_'
        else:
            print("Wrong scheme selected")
        return figtitle,figname
    
    def solution(self,CFL0,deltaX0,ui0,tf0,a0,N0,NP0,scheme):
        """
        Function: General function to solve advection equation using any type of scheme

        Args:
            CFL0: CFL number
            deltaX0: Grid Size
            ui0: Initial conditions
            tf0: Final time
            a0: Speed of wave
            N0: Number of Mesh points
            NP0: pllotting location in terms of number of time steps
            Scheme: Scheme
        """
        deltaT=CFL0*deltaX0/a0; t0=0; tplot=0;
        usol=copy.deepcopy(ui0); usoln=copy.deepcopy(usol); ufinal=np.zeros((N0,1));
        ufinal[1:N0-1,:]=usoln; ufinal[0,0]=0; ufinal[N0-1,0]=0;
        figtitle,figname=self.assemblematrix(scheme,CFL0,N0)
        umax=0;  umin=0;

        fig=plt.figure();
        plt.title(figtitle+'for  CFL='+str(CFL0));
        plt.xlabel('x')
        plt.ylabel('u')

        while t0-tf0<1e-6:
            x=np.linspace(0,1,N0)
            if abs(t0-tplot)<1e-5:
                tplot=t0+NP0*deltaT                
                max_val=np.max(ufinal)
                I=np.argmax(ufinal[:,0])

                plt.plot(x,ufinal,color='b',LineWidth=2)
                plt.text(x[I],max_val+0.01,'t='+str(format(t0, '.2f')))

                umax=max(umax,np.max(ufinal[:,0]))
                umin=min(umin,np.min(ufinal[:,0]))
            
            usoln=self.solnscheme.evaluate(N0,usol,usoln)

            usol=copy.deepcopy(usoln)
            ufinal[1:N0-1,0]=usoln[:,0]
            t0=t0+deltaT
            plt.xlim([0,1])
            plt.ylim([umin-0.1,umax+0.15])

        plt.grid(True,lw = 0.5)
        print("Close the figure with title",figtitle,"for next solution")
        plt.show()
        fig.savefig(figname+'for_CFL_'+str(CFL0)+'.png');
        plt.close()

if __name__ == '__main__':
    """
    Solution usin gall schemes
    """

    a0=1; #Speed of Wave
    L0=1; #Length of the Interval
    deltaX=0.005; #Grid Size
    CFL=[0.4, 1.0, 1.3]; #Required CFL Numbers to run
    N=int(L0/deltaX+1) #Number of grid points in the interval
    xf=0.75; #required position of the discontinuity
    NP=40; #Number of timesteps to plot the solution

    #Initial Conditions
    ui=np.zeros((N-2,1))
    for i in range(0,N-2):
        x=(i-1)*deltaX;
        if x>0.2 and x<=0.3:
            ui[i,0]=1   
    tf=(xf-0.3)/a0;

    lwe_obj=LinearWaveEqn()
    
    #Calling First Order Backward Difference scheme
    lwe_obj.solution(CFL[0],deltaX,ui,0.45,a0,N,round(tf/5/(CFL[0]*deltaX/a0)),scheme="FOBD_explicit");
    lwe_obj.solution(CFL[1],deltaX,ui,0.45,a0,N,round(tf/2/(CFL[1]*deltaX/a0)),scheme="FOBD_explicit")
    lwe_obj.solution(CFL[2],deltaX,ui,0.2,a0,N,round(0.1/1/(CFL[1]*deltaX/a0)),scheme="FOBD_explicit")

    #Calling Second Order Backward Difference scheme
    lwe_obj.solution(CFL[0],deltaX,ui,0.1,a0,N,round(tf/10/(CFL[0]*deltaX/a0)),scheme="SOBD_explicit")
    lwe_obj.solution(CFL[1],deltaX,ui,0.08,a0,N,round(tf/10/(CFL[1]*deltaX/a0)),scheme="SOBD_explicit")
    lwe_obj.solution(CFL[2],deltaX,ui,0.08,a0,N,round(tf/10/(CFL[2]*deltaX/a0)),scheme="SOBD_explicit")

    #Calling Crank-Nicolson Implicit scheme
    lwe_obj.solution(CFL[0],deltaX,ui,0.2,a0,N,round(tf/5/(CFL[0]*deltaX/a0)),scheme="CN_implicit")  
    lwe_obj.solution(CFL[1],deltaX,ui,0.2,a0,N,round(tf/5/(CFL[1]*deltaX/a0)),scheme="CN_implicit")
    lwe_obj.solution(CFL[2],deltaX,ui,0.2,a0,N,round(tf/5/(CFL[2]*deltaX/a0)),scheme="CN_implicit")
    
    #Calling Lax-Wendroff Scheme
    lwe_obj.solution(CFL[0],deltaX,ui,0.45,a0,N,round(0.45/5/(CFL[0]*deltaX/a0)),scheme="LW_explicit");  
    lwe_obj.solution(CFL[1],deltaX,ui,0.45,a0,N,round(0.45/2/(CFL[1]*deltaX/a0)),scheme="LW_explicit"); 
    lwe_obj.solution(CFL[2],deltaX,ui,0.1,a0,N,round(0.1/1/(CFL[2]*deltaX/a0)),scheme="LW_explicit");

    #Calling Lax Scheme
    lwe_obj.solution(CFL[0],deltaX,ui,0.45,a0,N,round(tf/5/(CFL[0]*deltaX/a0)),scheme="LAX_explicit")  
    lwe_obj.solution(CFL[1],deltaX,ui,0.45,a0,N,round(tf/2/(CFL[1]*deltaX/a0)),scheme="LAX_explicit") 
    lwe_obj.solution(CFL[2],deltaX,ui,0.2,a0,N,round(tf/5/(CFL[2]*deltaX/a0)),scheme="LAX_explicit")  

    #Calling First order implicit scheme
    lwe_obj.solution(CFL[0],deltaX,ui,0.45,a0,N,round(0.45/5/(CFL[0]*deltaX/a0)),scheme="FOBD_implicit")
    lwe_obj.solution(CFL[1],deltaX,ui,0.45,a0,N,round(0.45/5/(CFL[1]*deltaX/a0)),scheme="FOBD_implicit")  
    lwe_obj.solution(CFL[2],deltaX,ui,0.45,a0,N,round(0.45/5/(CFL[2]*deltaX/a0)),scheme="FOBD_implicit")
