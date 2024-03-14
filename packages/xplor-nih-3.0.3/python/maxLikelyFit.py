class MaxLikelyFit:
   """                                                                       
   A class which implements Maximum Likelihood Superpositioning of
   protein structures.
   Written by
           Robin Augustine Thottugal
           08/08/2007

  based on: 
    Theobald D.L. and Wuttke D.S., Empirical Bayes hierarchical
    models for regularizing maximum likelihood estimation in the matrix
    Gaussian Procrustes problem, Proc. of the National Academy of
    Sciences 103, 18521-18527 (2006a)

    Theobald D.L. and Wuttke D.S., THESEUS: Maximum likelihood
    superpositioning and analysis of macromolecular structures,
    Bioinformatics 22, 2171-2172 (2006)
     
   
   """
   
   def __init__(s,
                atmSelect,
                verbose, 
                mSelect=1):
      """
      atmSelect   - specify which atoms will be used in the fit
      mSelect     - specify method of regularizing the 
                    inverse of the covariance matrix -
                    if 0 use perturbative approach - add a small value to
                    the diagonal. If 1 use diagonal elements as ``eigenvalues,'' 
                    and follow the scheme of ref. FIX ME.
      verbose     - if 1, print out intermediate values of log likelihood and
                    RMSD.
      """
      s.atmSelect=str(atmSelect)
      s.mSelect=mSelect
      s.verbose=verbose
   

   def vecToMat(s,Vector):
      """
      Fuction which takes a cdsVector and
      convert it into a matrix 
      """
      from cdsMatrix import CDSMatrix_double as Mat
      M=Mat(len(Vector),3,0)
      atm_num=-1
      for elem in Vector:
         atm_num+=1   
         M[atm_num,0]=elem[0]
         M[atm_num,1]=elem[1]
         M[atm_num,2]=elem[2]
      return M

   def MatToVec(s,Structure):
      """
      Function which takes a cdsMatrix and
      convert into a cdsVector
      """
      from cdsVector import CDSVector_Vec3 as Vec
      FinalMeanStruct=Vec(Structure.rows())
      from vec3 import Vec3
      for i in range(Structure.rows()):
         x=Structure[i,0]
         y=Structure[i,1]
         z=Structure[i,2]
         Coord=Vec3(x,y,z)
         FinalMeanStruct[i]=Coord
      return FinalMeanStruct 


   def atmToMat(s,atm_list):
      """
      Function which converts a list of atom object to a cdsMatrix
      
      """
      from atomSel import AtomSel
      from cdsMatrix import CDSMatrix_double as Mat
      M=Mat(len(atm_list),3,0)
      atm_num=0
      for elem in atm_list:
         M[atm_num,0]=elem.pos()[0]
         M[atm_num,1]=elem.pos()[1]
         M[atm_num,2]=elem.pos()[2]
         atm_num+=1 
      return M   


   def cnorm(s,M,pM):
      """
      Calculates the RMSD between two matrix
      """
      from math import sqrt
      nm=0
      from cdsVector import CDSVector_double as V
      from cdsVector import norm
      Vect_1=[V]*M.rows()
      Vect_2=[V]*M.rows()
      for i in range(M.rows()):
         Vect_1[i]=V(M.cols())
         Vect_2[i]=V(M.cols()) 
         for j in range(M.cols()):
            Vect_1[i][j]=M[i,j]
            Vect_2[i][j]=pM[i,j]
         dis=norm(Vect_1[i]-Vect_2[i])
         dis=(dis*dis)
         nm+=dis    
      return (nm/M.rows()) 

   def newton(s,x,Mean_inv,inv_egnval):
      """
      Newton Rapson Method for root finding.This function is
      used to find gamma
      """
      from math import log
      from moremath import digamma,trigamma
      Th_hold=1e-6    
      log_v=0
      for i in range(len(inv_egnval)):
         log_v+=log(inv_egnval[i])
      log_val=log_v/len(inv_egnval)
      constant=-log(Mean_inv)+log_val
      count = 0
      while 1:
         f=log(x)-digamma(x)+constant
         fd=(1.0/x)-trigamma(x)
         dx =f/float(fd)
         xn=x-dx
         error=((xn-x)/xn)*100
         if (abs(error) < Th_hold):
            return (x,error,f)
         if(xn<1e-6):
            x=1e-6
            continue
         x=x-dx   
         count = count + 1
         #print "newton:", count, abs(error)


   def translate(s,CovMat,X,NStruct):
      """
      Function to find the ML estimate of optimal translation
      """
      from cdsMatrix import CDSMatrix_double as Mat    
      from cdsMatrix import inverse,transpose
      One_k=Mat(CovMat.rows(),1,1)
      transVec=list(range(NStruct))
      for i in range(NStruct):
         transVec[i]=Mat(1,3,0)
      X_cap=list(range(NStruct))
      den=transpose(One_k)*inverse(CovMat)*One_k
      for i in range(NStruct):
         num=transpose(One_k)*inverse(CovMat)*X[i]
         transVec[i]=(-1.0/den[0,0])*num
         X_cap[i]=X[i]+(One_k*transVec[i])
         pass
      return (X_cap,transVec)

   def rotate(s,MStruct,CovMat,X_cap,NStruct):
      """
      Function to find the ML estimate of optimal Rotation 
      """
      from cdsMatrix import CDSMatrix_double as Mat
      from cdsMatrix import inverse,transpose,svd,det
      RotMat=list(range(NStruct))
      for i in range(NStruct):
         RotMat[i]=Mat(3,3,0)
      X_Rot=list(range(NStruct)) 
      prev_MStruct=Mat(MStruct.rows(),MStruct.cols())
      prev_MStruct=MStruct
      while(1):
         Avg=Mat(MStruct.rows(),MStruct.cols(),0)
         for i in range (NStruct):
            SVD_Mat=svd(transpose(MStruct)*inverse(CovMat)*X_cap[i])
            P=Mat(len(SVD_Mat.sigma),len(SVD_Mat.sigma),0)
            P.setDiag(1)
            val=det(SVD_Mat.u)*det(SVD_Mat.vT)
            if (val.real<0):
               P[2,2]=-1
            RotMat[i]=transpose(SVD_Mat.vT)*P*transpose(SVD_Mat.u)
            if(det(RotMat[i]).real<0):
               print('Warning:Determinant of Rotation Matrix for Structure[%d]is %f'%(i,det(RotMat[i]).real))
            X_Rot[i]=X_cap[i]*RotMat[i]
            Avg+=X_Rot[i]
            pass   
         MStruct=(1.0/NStruct)*Avg
         d=s.cnorm(MStruct,prev_MStruct)
         if(d<0.001):
            break
         prev_MStruct=MStruct  
         pass
      return(RotMat,X_Rot,MStruct)           

   def estimateEgien(s,Segnval,egnval,NStruct):
      """
      Function which implements the Expectation Maximization algorithm to
      fit all the eigen value according to an inverse gamma distribution
      """
      invegnval=list(range(len(egnval)))
      prev_gamma=0
      prev_alpha=0
      inum=1
      while(1):
         avg=0.0
         var=0.0
         for i in range(len(egnval)):
            invegnval[i]=1.0/egnval[i]
            avg+=invegnval[i]
            var+=(invegnval[i]*invegnval[i])
            pass
         Mean_inv=avg/len(egnval)
         Variance=(float(var)/len(egnval))-(Mean_inv*Mean_inv)
         gamma=(Mean_inv*Mean_inv)/Variance
         alpha=Mean_inv/Variance
         if(inum==1):
            for i in range(len(egnval)):
               num=(3.0*NStruct)*(((2.0*alpha)/(3.0*NStruct))+Segnval[i])
               den=(3.0*NStruct)+(2.0*(gamma+1))
               egnval[i]=float(num)/den
               pass
            inum=0
            continue
         result=s.newton(gamma,Mean_inv,invegnval)
         new_gamma=result[0]
         new_alpha=float(new_gamma)/Mean_inv
         if (abs(prev_gamma-new_gamma)<1e-8):
            return(new_gamma,new_alpha,result[2])
         for i in range(len(egnval)):
            num=(3.0*NStruct)*(((2.0*new_alpha)/(3.0*NStruct))+Segnval[i])
            den=(3.0*NStruct)+(2.0*(new_gamma+1))
            egnval[i]=float(num)/den
            pass
         prev_gamma=new_gamma                          
              

   def loglikelihood(s,CovMat,X,MStruct,RotMat,transVec,NStruct,gamma,alpha):
      """
      Function which evaluvates the log likelihood function 
      """
      from cdsMatrix import CDSMatrix_double as Mat
      from cdsMatrix import eigen,trace,det,inverse,transpose
      from math import log
      from moremath import gamma as Gamma
      One_k=Mat(CovMat.rows(),1,1)
      K=CovMat.rows()
      norm=0
      list_egn=eigen(CovMat)
      deter=0
      Y=list(range(NStruct))
      deter=0
      for i in range(CovMat.rows()):
         deter+=log(CovMat[i,i])
      for i in range(NStruct):
         Y[i]=(X[i]+(One_k*transVec[i]))*RotMat[i]
         norm+=trace(transpose(Y[i]-MStruct)*inverse(CovMat)*(Y[i]-MStruct))
      term_1=norm/2
      term_2=(3*NStruct*K*log(2*3.14))/2
      term_3=3*NStruct*deter/2
      term_4=(1+gamma)*deter
      term_5=alpha*trace(inverse(CovMat))
      term_6=K*gamma*log(alpha)
      term_7=K*log(Gamma(gamma))
      #Based on the method selected the loglikehood is evaluvated
      if(s.mSelect==1):
         ls=-term_1-term_2-term_3-term_4-term_5+term_6-term_7
      else:   
         ls=-term_1-term_3
      return ls
      pass
 
   def calcVariance(s,fitted):
      """
      Function which calculates the variance of all atoms
      which goes into the Bfactor column of the average
      structure.
      """
      from cdsMatrix import CDSMatrix_double as Mat
      from math import sqrt
      from cdsVector import norm
      import xplor
      import simulation
      var=[]
      sim=xplor.simulation
      NStruct = len(fitted)
      sim.setAtomPosArr(fitted[0])
      MStructure=sim.atomPosArr()
      for i in range(1,NStruct):
         sim.setAtomPosArr(fitted[i])
         MStructure+=sim.atomPosArr()
         pass
      MStructure/=NStruct
      MeanStr=s.vecToMat(MStructure)
      for i in range(len(fitted[0])):
         s=0
         for j in range(NStruct):
            x=(fitted[j][i][0]-MeanStr[i,0])
            y=(fitted[j][i][1]-MeanStr[i,1])
            z=(fitted[j][i][2]-MeanStr[i,2])
            s+=(x*x)+(y*y)+(z*z)
         s/=(NStruct)
         s=sqrt(s)
         var.append(s)
      return(var)
#   def calcVariance(s,Y,transVec,RotMat,NStruct):
#      """
#      Function which calculates the variance of all atoms
#      which goes into the Bfactor column of the average
#      structure.
#      """
#      from cdsMatrix import CDSMatrix_double as Mat
#      from math import sqrt
#      from cdsVector import norm
#      import xplor
#      import simulation
#      var=[]
#      sim=xplor.simulation
#      Z=range(NStruct)
#      One_k=Mat(Y[0].rows(),1,1)
#      Z[0]=(Y[0]+(One_k*transVec[0]))*RotMat[0]
#      sim.setAtomPosArr(s.MatToVec((Y[0]+(One_k*transVec[0]))*RotMat[0]))
#      MStructure=sim.atomPosArr()
#      for i in range(1,NStruct):
#         Z[i]=(Y[i]+(One_k*transVec[i]))*RotMat[i]
#         sim.setAtomPosArr(s.MatToVec((Y[i]+(One_k*transVec[i]))*RotMat[i]))
#         MStructure+=sim.atomPosArr()
#         pass
#      MStructure/=NStruct
#      MeanStr=s.vecToMat(MStructure)
#      for i in range(Y[0].rows()):
#         s=0
#         for j in range(NStruct):
#            x=(Z[j][i,0]-MeanStr[i,0])
#            y=(Z[j][i,1]-MeanStr[i,1])
#            z=(Z[j][i,2]-MeanStr[i,2])
#            s+=(x*x)+(y*y)+(z*z)
#         s/=(NStruct)
#         s=sqrt(s)
#         var.append(s)
#      return(var)
#
   
   def RMSD(s,X_Rot,MStruct,NStruct,Natoms):
      """
      Calculates the variance of selected atoms
      """
      from math import sqrt
      from vec3 import Vec3, norm
      
      s.RMSD_atm=list(range(Natoms))
      for i in range(Natoms):
         sum=0
         for j in range(NStruct):
            d=Vec3(X_Rot[j][i,0]-MStruct[i,0],
                   X_Rot[j][i,1]-MStruct[i,1],
                   X_Rot[j][i,2]-MStruct[i,2])
            sum+=norm(d)**2
            pass
         s.RMSD_atm[i]=sqrt(sum/NStruct)
         pass
      return
   
   def NormalRMSD(s,X_Rot,MStruct,NStruct,Natoms):
      """
      Calculates the Average RMSD of N Structure
      """
      from math import sqrt
      from vec3 import Vec3, norm
      nRMSD=0
      for i in range(NStruct):
         sum=0
         for j in range(Natoms):
            d=Vec3(X_Rot[i][j,0]-MStruct[j,0],
                   X_Rot[i][j,1]-MStruct[j,1],
                   X_Rot[i][j,2]-MStruct[j,2])
            sum+=norm(d)**2
            pass
         nRMSD+=sqrt(sum/Natoms)
         pass
      return(nRMSD/NStruct)
                   
   def fit(s,coords):
      """
      perform the maximum likelihood fitting procedure.
      """

      s.coords0 = coords
      def write(str):
         import sys
         sys.stdout.write(str)
         sys.stdout.flush()
         pass

      import simulation
      import xplor
      from atomSel import AtomSel
      from atomSelAction import Fit
      from math import sqrt
      from cdsMatrix import eigen,inverse,trace
      from cdsMatrix import CDSMatrix_double as Mat
      from cdsMatrix import transpose,det
       
      NStruct=len(coords)
      sim=xplor.simulation
      X=[]
#      s.Y=[]
      fitTo=sim.atomPosArr()
      MeanStruct=sim.atomPosArr()
      MeanStruct.set( [0,0,0] )
      fitsel=s.atmSelect
      for coord in coords:
#         s.Y.append( s.vecToMat(coord) )
         xplor.simulation.setAtomPosArr( coord )
         select=AtomSel(s.atmSelect) 
         X.append( s.atmToMat(select) )
         AtomSel("known").apply(Fit(fitTo,fitsel))
         MeanStruct+=sim.atomPosArr()
         pass
 
      MeanStruct/= NStruct
      s.SelectedMStruct=MeanStruct
      sim.setAtomPosArr(MeanStruct)

      #MStruct conatins only the selected atoms
      MStruct=s.atmToMat(select)
      #Intial Guess for the Covariance Matrix
      
      Natoms=len(select)
      s.SelectedAtoms=list(select)
      CovMat=Mat(Natoms,Natoms,0)
      CovMat.setDiag(1)
      s.i_number=0
      pMStruct=Mat(Natoms,3,0)
      while(1):
         s.i_number+=1
         #Calculating the Optimal translation'
         (X_cap,s.transVec)=s.translate(CovMat,X,NStruct)
         #'Calculating the Optimal Rotation'
         (s.RotMat,X_Rot,MStruct)=s.rotate(MStruct,CovMat,X_cap,NStruct)
         #The loop exits when the RMSD between the previous Mean structure and
         #current Mean Structure is small
         conv=abs(s.cnorm(pMStruct,MStruct)) 
         if (conv<1e-7):
            s.RMSD(X_Rot,MStruct,NStruct,Natoms)
            s.fittedStructs=X_Rot
            s.meanStruct = MStruct
            break
         pMStruct=MStruct
         #Calculating the Sample Co-variance matrix
         SCovMat=Mat(Natoms,Natoms,0)
         for i in range(NStruct):
            SCovMat+=(X_Rot[i]-MStruct)*transpose(X_Rot[i]-MStruct)
            pass
        
         SCovMat=(1.0/(3*NStruct))*SCovMat
         #Checking the Covariance matrix to make sure that not all diagonal elements
         #are close to zero.If all diagonal elements are zero the program exit 
         count=0
         for i in range(Natoms):
            for j in range(Natoms):
               if(i==j):
                  if((abs(SCovMat[i,i]))<1e-8):
                     count+=1
                  continue
               SCovMat[i,j]=0
         if(count==Natoms):
            import sys
            print('All Diagonal elements of the Covariance Matrix is Zero')
            sys.exit(0)
      
         Segnval=[]
         L=[]
      
         
         #Eigen Values Based on Eigen Decompostion
         #Li=[]
         #L=eigen(SCovMat)
         #for i in range(len(L)):
         #   Li.append(L[i].value().real)
         #   pass
         #Li.sort()
         #Th_hold=(Li[2]*10)
         #print Th_hold
         #for i in range(len(Li)):
         #   if(Li[i]>Th_hold):
         #       Segnval.append(Li[i])

                  
         if (s.mSelect==1):
            #Eigen Values Based on diagonal Elements
            for i in range(Natoms):
               if((abs(SCovMat[i,i])-0)<1e-8):
                  continue
               L.append(SCovMat[i,i])
               pass
            L.sort() 
            Segnval=list(L)
            egnval=list(Segnval)
            (s.gamma,s.alpha,fval)=s.estimateEgien(Segnval,egnval,NStruct)
            Identity_Mat=Mat(Natoms,Natoms,0)
            Identity_Mat.setDiag(1)
            term_1=((2*s.alpha)/(3*NStruct))*Identity_Mat
            term_1=term_1+SCovMat
            term_2=(3*NStruct)/((3*NStruct)+(2*(s.gamma+1)))
            SCovMat=term_2*term_1
         else:
            #adding a small value to the diagonal of Covarinace matrix to make
            #sure that an inverse exist for the matrix
            s.gamma=0.01
            s.alpha=0.01
            Identity_Mat=Mat(Natoms,Natoms,0)
            Identity_Mat.setDiag(1e-4)
            SCovMat+=Identity_Mat
            pass
            
         CovMat=SCovMat
         s.RMSD_ml=sqrt(Natoms/(trace(inverse(CovMat))))
         s.LsSigma=sqrt((trace(CovMat))/Natoms)
         s.ls=s.loglikelihood(CovMat,X,MStruct,s.RotMat,s.transVec,
                              NStruct,s.gamma,s.alpha)

         s.nRMSD=s.NormalRMSD(X_Rot,MStruct,NStruct,Natoms)
         if s.verbose==0:
            write('.')
         
         if s.verbose==1:
            print('Iteration Number:%d'%s.i_number)
            print('  RMSD_ml: %.3f'%s.RMSD_ml)
            print('  Least Square Sigma: %.3f'%s.LsSigma)
            print('  Log Likelihood: %.2f'%(s.loglikelihood(CovMat,X,
                                                            MStruct,s.RotMat,
                                                            s.transVec,NStruct,
                                                            s.gamma,s.alpha)))     

   def WriteRotTransMat(s,NStruct):
       'Writing the Rotation and Translation Matrix'
       f=open('RotTrans_Mat.txt','w')
       for i in range(NStruct):
          f.write('\n Translation and Rotation Matrix for Structure:%d'%i)
          f.write('\n Translation Matrix')
          f.write('\n \t %f \t %f \t %f'%(s.transVec[i][0,0],
                                          s.transVec[i][0,1],
                                          s.transVec[i][0,2]))
          f.write('\n Rotation Matrix')
          for j in range(3):
             f.write('\n')
             for k in range(3):
                f.write('\t %f'%s.RotMat[i][j,k])
                pass
             pass
          pass
       return
 
   def PDBWrite(s,structs,wStruct):
      """Writing the Co-ords into PDB files

      if averageFilename is specified,
      average coordinates will also be written to the specified filename.
      """
      import xplor
      from atomSel import AtomSel
      from atomSelAction import Rotate, Translate
      import simulation
      from pdbTool import PDBTool
      from cdsMatrix import CDSMatrix_double as Mat
      from vec3 import Vec3
      from mat3 import Mat3, transpose
      sim=xplor.simulation
#      One_k=Mat(s.Y[0].rows(),1,1)
      fitted=[]
      sim.setAtomPosArr(s.coords0[0])
      AtomSel("all").apply( Translate(Vec3(s.transVec[0][0,0],
                                           s.transVec[0][0,1],
                                           s.transVec[0][0,2])) )
      rot=Mat3(s.RotMat[0][0,0],s.RotMat[0][0,1],s.RotMat[0][0,2],
                    s.RotMat[0][1,0],s.RotMat[0][1,1],s.RotMat[0][1,2],
                    s.RotMat[0][2,0],s.RotMat[0][2,1],s.RotMat[0][2,2])
      AtomSel("all").apply(Rotate(transpose(rot)))
      fitted.append( sim.atomPosArr())
#s.MatToVec((s.Y[0]+(One_k*s.transVec[0]))*
#                                   s.RotMat[0] ))
      MStructure=sim.atomPosArr()
      if(wStruct==1):
         p_obj=PDBTool("mle_"+structs[0])
         p_obj.write()
      for i in range(1,len(structs)):
#         sim.setAtomPosArr(s.MatToVec((s.Y[i]+(One_k*s.transVec[i]))*
#                                      s.RotMat[i]))
         sim.setAtomPosArr(s.coords0[i]) 
         rot = Mat3(s.RotMat[i][0,0],s.RotMat[i][0,1],s.RotMat[i][0,2],
                    s.RotMat[i][1,0],s.RotMat[i][1,1],s.RotMat[i][1,2],
                    s.RotMat[i][2,0],s.RotMat[i][2,1],s.RotMat[i][2,2])  
         AtomSel("all").apply( Translate(Vec3(s.transVec[i][0,0],
                                           s.transVec[i][0,1],
                                           s.transVec[i][0,2])) )
         AtomSel("all").apply(Rotate(transpose(rot)))
         fitted.append( sim.atomPosArr())
         MStructure+=sim.atomPosArr()
         if(wStruct==1):
            p_obj=PDBTool('mle_'+structs[i])
            p_obj.write()
         pass

#     var=s.calcVariance(s.Y,s.transVec,s.RotMat,len(structs))
        
      var=s.calcVariance(fitted)
      MStructure/=len(structs)
      
      sim.setAtomPosArr(s.SelectedMStruct)
      s.variance=[]
              
      for i in range(sim.numAtoms()):
          s.variance.append((sim.atomByID(i),var[i])) 
      if(wStruct==1):  
         sim.setAtomPosArr(MStructure)
         p_obj=PDBTool("MLE-Average.pdb")
         for i in range(sim.numAtoms()):
            p_obj.setAux2(sim.atomByID(i),var[i])
            pass
         p_obj.write()
      return
   def fittedCoords(s):
      """
      return an array of all atomic coordinates with the most recent
      fit translate+rotation applied
      """
      from atomSel import AtomSel
      from atomSelAction import Translate, Rotate
      from vec3 import Vec3
      from mat3 import Mat3, transpose
      import xplor
      fitted=[]
      sim=xplor.simulation
      for i in range(len(s.coords0)):
         sim.setAtomPosArr( s.coords0[i] )
         AtomSel("all").apply( Translate(Vec3(s.transVec[i][0,0],
                                              s.transVec[i][0,1],
                                              s.transVec[i][0,2])))
         rot = Mat3(s.RotMat[i][0,0],s.RotMat[i][0,1],s.RotMat[i][0,2],
                    s.RotMat[i][1,0],s.RotMat[i][1,1],s.RotMat[i][1,2],
                    s.RotMat[i][2,0],s.RotMat[i][2,1],s.RotMat[i][2,2])
         AtomSel("all").apply(Rotate(transpose(rot)))
         fitted.append( sim.atomPosArr() )
         pass

      return fitted
      
   def orderedResidues(s,threshold,
                       selection=None):
      """
      return a tuple describing the range of ordered residues.
      The first element is the total number of ordered residues, and the
      second element is a list of tuples corresponding to ranges of
      orderedResidues.

      A residue is defined as ordered if the average positional RMSD of
      selected atoms of that residue are less than threshold.

      If selection is not specified, self.atmSelect will be used.
      """

      if not selection:
         selection = s.atmSelect
         pass

      if type(selection)==type("string"):
         import atomSel
         selection=atomSel.AtomSel(selection)
         pass

      fitted=s.fittedCoords()
         
      from cdsVector import CDSVector_Vec3, CDSVector_double, vec_norm, sqrt
      mean=CDSVector_Vec3(fitted[0])
      for coord in fitted[1:]:
         mean += coord
      mean /= len(fitted)

      rmsd=CDSVector_double(len(mean),0.)
      for coord in fitted:
         rmsd += vec_norm(coord-mean)**2
         pass
      rmsd /= len(fitted)
      rmsd = sqrt(rmsd)
      

      residues={}
      cnt=0
      for atom in selection:
         resid = atom.residueNum()
         if resid not in residues:
            residues[resid]=[0,0.] #num atoms, cumulative rmsd
            pass
         residues[resid][0] += 1
#         print atom.string(),
#         print "%6.2f %6.2f" %(rmsd[ atom.index() ],
#                               s.RMSD_atm[cnt]),
#         for i in range(len(fitted)):
#            print fitted[i][atom.index()],
#            print "(%6.2f %6.2f %6.2f)" % (s.fittedStructs[i][cnt,0],
#                                           s.fittedStructs[i][cnt,1],
#                                           s.fittedStructs[i][cnt,2]),
#            pass
#         print
         residues[resid][1] += rmsd[ atom.index() ]
         cnt+=1
         pass

      ordered=[]
      num=0
      for resid in list(residues.keys()):
         aveRMSD=residues[resid][1] / residues[resid][0]
         if aveRMSD<threshold:
            ordered.append(resid)
            num += 1
            pass
         pass

      resids=ordered
      resids.sort()

      seq=[]
      cnt=0
      while cnt<len(resids):
         first=resids[cnt]
         while (cnt+1<len(resids)        and
                resids[cnt]+1==resids[cnt+1]):
            cnt += 1
            pass
         last=resids[cnt]
         cnt += 1
         seq.append( (first,last) )
         pass

      return (num,seq)
      
      
         
