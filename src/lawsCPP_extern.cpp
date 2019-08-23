#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Eigen;
using namespace Rcpp;
using namespace std;


// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
Eigen::MatrixXd AAt(const Eigen::MatrixXd& A) {
    using   Eigen::LLT;
    using   Eigen::Lower;
    using   Eigen::Map;
    using   Eigen::MatrixXd;
    using   Eigen::MatrixXi;
    using   Eigen::Upper;
    using   Eigen::VectorXd;
    using   Eigen::VectorXi;
    //typedef Map<MatrixXd>  MapMatd;
    //typedef Map<MatrixXi>  MapMati;
    //typedef Map<Eigen::VectorXd>  MapVecd;
    //typedef Map<VectorXi>  MapVeci;
    
    int    n(A.cols());
    return   MatrixXd(n,n).setZero().selfadjointView<Lower>()
                          .rankUpdate(A);
}

// [[Rcpp::export]]
Eigen::MatrixXd AtA(const Eigen::MatrixXd& A) {
    using   Eigen::LLT;
    using   Eigen::Lower;
    using   Eigen::Map;
    using   Eigen::MatrixXd;
    using   Eigen::MatrixXi;
    using   Eigen::Upper;
    using   Eigen::VectorXd;
    using   Eigen::VectorXi;
    //typedef Map<MatrixXd>  MapMatd;
    //typedef Map<MatrixXi>  MapMati;
    //typedef Map<Eigen::VectorXd>  MapVecd;
    //typedef Map<VectorXi>  MapVeci;
    
    int    n(A.cols());
    return   MatrixXd(n,n).setZero().selfadjointView<Lower>()
                          .rankUpdate(A.adjoint());
}

// [[Rcpp::export]]
List lltLS(const    Eigen::Map<Eigen::VectorXd>&    y,
           const   Eigen::MatrixXd&   B,
                  const   double&     tau,
                  const   Eigen::Map<Eigen::VectorXd>&    lambdashort_glatt,
                  const   Eigen::Map<Eigen::VectorXd>&    lambdashort_orig,
                  Eigen::MatrixXd    DD,
                  const   Eigen::Map<Eigen::VectorXi>&    NB,
                  const   Eigen::Map<Eigen::VectorXi>&    glatterms) {
    
    using   Eigen::LLT;
    using   Eigen::Lower;
    using   Eigen::Map;
    using   Eigen::MatrixXd;
    using   Eigen::MatrixXi;
    using   Eigen::Upper;
    using   Eigen::VectorXd;
    using   Eigen::VectorXi;
    //typedef Map<MatrixXd>  MapMatd;
    //typedef Map<MatrixXi>  MapMati;
    //typedef Map<Eigen::VectorXd>  MapVecd;
    //typedef Map<VectorXi>  MapVeci;
    
    VectorXd lambdalong(NB.sum()+1);
    lambdalong.fill(1);
    int lli = 0;
    for(int li=0;li<B.cols()-NB.sum();li++){
        lambdalong[li] = 0;
        lli++;
    }
    VectorXd lambdashort = lambdashort_orig;
    
    for(int k=0; k < glatterms.size(); k++) {
        if(NB(glatterms(k)-1) != 0) {
            lambdashort(glatterms(k)-1) = lambdashort_glatt(k);
        } 
    }
    for(int nbi=0;nbi<NB.size();nbi++){
        for(int li=0;li<NB[nbi];li++){
            lambdalong[lli] = lambdashort[nbi];
            lli++;
        }
    }
    
    for(int i=0; i<DD.cols();i++)
        DD.col(i) *= sqrt(lambdalong[i]);
    
    MatrixXd X (B.rows()+DD.rows(),B.cols());
    X << B,DD;
    
    
    
    
    const int n(B.rows());
    
    LLT<MatrixXd> llt(AtA(X));
    
    VectorXd  betahat(llt.solve(X.adjoint() * y));
    
    VectorXd  fitted(X * betahat);
    
    VectorXd weight(X.rows());
    
    weight.fill(1);
    
    
    for(int i=0;i<n;i++)
        if(y[i] < fitted[i])
            weight[i] = 1 - tau;
        else
            weight[i] = tau;
        
        //return List::create(Named("weights") =weight);
        
        //double dw = 1;
        int iter = 1;
        VectorXd wold = weight;
        wold[1]++;
        while((wold - weight).norm() != 0 && iter < 50)
        {
            checkUserInterrupt();
            
            llt.compute(X.adjoint() * weight.asDiagonal() * X);
            betahat = llt.solve(X.adjoint() * weight.asDiagonal() * y);
            
            iter++;
            
            wold = weight;
            fitted = X*betahat;
            
            for(int i=0;i<n;i++)
                if(y[i] < fitted[i])
                    weight[i] = 1 - tau;
                else
                    weight[i] = tau;
                
        }
        
        return List::create(Named("a") =betahat,Named("weight") =weight.head(n),Named("fitted") = fitted.head(n),Named("diag.hat.ma") = (X*llt.solve(X.adjoint() * weight.asDiagonal())).diagonal().head(n));
}


// [[Rcpp::export]]
double smoothCPPP(const    Eigen::Map<Eigen::VectorXd>&    y,
                  const   Eigen::MatrixXd&   B,
                  const   double&     tau,
                  const   Eigen::Map<Eigen::VectorXd>&    lambdashort_glatt,
                  const   Eigen::Map<Eigen::VectorXd>&    lambdashort_orig,
                  Eigen::MatrixXd    DD,
                  const   Eigen::Map<Eigen::VectorXi>&    NB,
                  const   Eigen::Map<Eigen::VectorXi>&    glatterms, 
                  std::string smoothtype){
    using   Eigen::LLT;
    using   Eigen::Lower;
    using   Eigen::Map;
    using   Eigen::MatrixXd;
    using   Eigen::MatrixXi;
    using   Eigen::Upper;
    using   Eigen::VectorXd;
    using   Eigen::VectorXi;
    //typedef Map<MatrixXd>  MapMatd;
    //typedef Map<MatrixXi>  MapMati;
    //typedef Map<Eigen::VectorXd>  MapVecd;
    //typedef Map<VectorXi>  MapVeci;
    
    VectorXd lambdalong(NB.sum()+1);
    lambdalong.fill(1);
    int lli = 0;
    for(int li=0;li<B.cols()-NB.sum();li++){
        lambdalong[li] = 0;
        lli++;
    }
    VectorXd lambdashort = lambdashort_orig;
    
    for(int k=0; k < glatterms.size(); k++) {
        if(NB(glatterms(k)-1) != 0) {
            lambdashort(glatterms(k)-1) = lambdashort_glatt(k);
        } 
    }
    for(int nbi=0;nbi<NB.size();nbi++){
        for(int li=0;li<NB[nbi];li++){
            lambdalong[lli] = lambdashort[nbi];
            lli++;
        }
    }
    /*Rcout << lambdalong << std::endl;
     Rcout << lambdashort << std::endl;
     Rcout << std::endl;*/
    
    for(int i=0; i<DD.cols();i++)
        DD.col(i) *= sqrt(lambdalong[i]);
    
    MatrixXd X (B.rows()+DD.rows(),B.cols());
    X << B,DD;
    
    
    const int n(B.rows());
    
    LLT<MatrixXd> llt(AtA(X));
    
    VectorXd  betahat(llt.solve(X.adjoint() * y));
    
    VectorXd  fitted(X * betahat);
    
    VectorXd weight(X.rows());
    
    weight.fill(1);
    
    
    for(int i=0;i<n;i++)
        if(y[i] < fitted[i])
            weight[i] = 1 - tau;
        else
            weight[i] = tau;
        
        //return List::create(Named("weights") =weight);
        
        //double dw = 1;
        int iter = 1;
        VectorXd wold = weight;
        wold[1]++;
        while((wold - weight).norm() != 0 && iter < 50)
        {
            checkUserInterrupt();
            
            llt.compute(X.adjoint() * weight.asDiagonal() * X);
            betahat = llt.solve(X.adjoint() * weight.asDiagonal() * y);
            
            iter++;
            
            wold = weight;
            fitted = X*betahat;
            
            for(int i=0;i<n;i++)
                if(y[i] < fitted[i])
                    weight[i] = 1 - tau;
                else
                    weight[i] = tau;
                
        }
        
        
        double score = 0;
        /*double tempscore1 = ((weight.head(n).array())*((pow((((y.head(n) - fitted.head(n)).array())),2)).array())).sum();
        double tempscore2 = log((((weight.head(n).array())*((pow((((y.head(n) - fitted.head(n)).array())),2)).array())).sum())/n)*n;
        double tempscore3 = 2*((X*llt.solve(X.adjoint() * weight.asDiagonal())).diagonal().head(n).sum());
        Rcout << tempscore1 << tempscore2 << tempscore3 << std::endl;*/
        if(smoothtype == "aic") {
            score = log((((weight.head(n).array())*((pow((((y.head(n) - fitted.head(n)).array())),2)).array())).sum())/n)*n + 2*((X*llt.solve(X.adjoint() * weight.asDiagonal())).diagonal().head(n).sum());
        }
        if(smoothtype == "bic") {
            score = log((((weight.head(n).array())*((pow((((y.head(n) - fitted.head(n)).array())),2)).array())).sum())/n)*n + log(n)*((X*llt.solve(X.adjoint() * weight.asDiagonal())).diagonal().head(n).sum());
        }
        if(smoothtype == "ocv") {
            score = ( (((weight.head(n).array())*((pow((((y.head(n) - fitted.head(n)).array())),2)).array()))).array() / 
                          (pow(((1-((X*llt.solve(X.adjoint() * weight.asDiagonal())).diagonal().head(n)).array()).array()),2)).array() ).mean();
        }
        if(smoothtype == "gcv") {
            score = (((weight.head(n).array())*((pow((((y.head(n) - fitted.head(n)).array())),2)).array())).mean()) * 1/pow(1-(1+(X*llt.solve(X.adjoint() * weight.asDiagonal())).diagonal().head(n).sum())/n,2);
        }
        
        if(!(score != score) && !(score==score))
        {
            score = 10^50;
        }
        
        return (score);
}


// [[Rcpp::export]]
List schallCPPfun(const     Eigen::Map<Eigen::VectorXi>&    glatterms, 
                  const   Eigen::Map<Eigen::VectorXd>&    y,
                  const   Eigen::MatrixXd&                B,
                  const   double&                         tau,
                  const   Eigen::Map<Eigen::VectorXd>&    lambdashort_in,
                  Eigen::MatrixXd                 DD,
                  const   Eigen::Map<Eigen::VectorXi>&    NB, 
                  const   Eigen::Map<Eigen::VectorXi>&    NBP, 
                  Eigen::Map<Eigen::VectorXi>     NBPC, 
                  bool                            center){
    
    using   Eigen::LLT;
    using   Eigen::Lower;
    using   Eigen::Map;
    using   Eigen::MatrixXd;
    using   Eigen::MatrixXi;
    using   Eigen::Upper;
    using   Eigen::VectorXd;
    using   Eigen::VectorXi;
    //typedef Map<MatrixXd>  MapMatd;
    //typedef Map<MatrixXi>  MapMati;
    //typedef Map<Eigen::VectorXd>  MapVecd;
    //typedef Map<VectorXi>  MapVeci;
    
    // NB   is the number of columns,coefficients, etc. per covariate
    // NBP  is the number of penalized coeffients per covariate, so NB-1 for 
    //      P-splines, NB for GMRF, and NB-3 for 2dsplines
    // NBPC is the number of effects, before the current covariates penalized 
    //      part starts (includes the position of the unpenalized part of the 
    //      current covariate) (only used in schall algorithm) 
    //      (position where to start block for penalized segment of current covariate (in R, in c++ -1))
    
    VectorXd lambdashort = lambdashort_in;
    int m = y.size();
    
    double dc = 1;
    double dw = 1;
    int it = 1;
    int iter = 1;
    
    VectorXd lambdalong(NB.sum()+1);
    lambdalong.fill(1);
    
    
    // # Start Initialize LAWS for given lambda
    int lli = 0;
    
    MatrixXd X (B.rows()+DD.rows(),B.cols());
    
    MatrixXd DDorig = DD;
    const int n(B.rows());
    
    
    VectorXd weight(X.rows());
    
    weight.fill(1);
    VectorXd w0 = weight;
    
    
    w0 = weight;
    DD = DDorig;
    
    lli = 0;
    for(int li=0;li<B.cols()-NB.sum();li++){
        lambdalong[li] = 0;
        lli++;
    }
    
    for(int nbi=0;nbi<NB.size();nbi++){
        for(int li=0;li<NB(nbi);li++){
            lambdalong[lli] = lambdashort[nbi];
            lli++;
        }
    }
    
    for(int i=0; i<DD.cols();i++)
        DD.col(i) *= sqrt(lambdalong[i]);
    
    X << B,DD;
    LLT<MatrixXd> llt(AtA(X));
    
    VectorXd betahat = (llt.solve(X.adjoint() * y));
    VectorXd fitted = (X * betahat);
    weight.fill(1);
    for(int i=0;i<n;i++)
        if(y[i] < fitted[i])
            weight[i] = 1 - tau;
        else
            weight[i] = tau;
        
        
        while((dc >= 0.01 || dw != 0) && it < 100) {
            checkUserInterrupt();
            
            // # Start Fit LAWS for given lambda
            if(it > 1) {
                w0 = weight;
                DD = DDorig;
                lli = 0;
                for(int li=0;li<B.cols()-NB.sum();li++){
                    lambdalong[li] = 0;
                    lli++;
                }
                
                for(int nbi=0;nbi<NB.size();nbi++){
                    for(int li=0;li<NB(nbi);li++){
                        lambdalong[lli] = lambdashort[nbi];
                        lli++;
                    }
                }
                
                for(int i=0; i<DD.cols();i++)
                    DD.col(i) *= sqrt(lambdalong[i]);
                
                X << B,DD;
                LLT<MatrixXd> llt(AtA(X));
                betahat = (llt.solve(X.adjoint() * y));
                
                fitted = (X * betahat);
                weight.fill(1);
                for(int i=0;i<n;i++)
                    if(y[i] < fitted[i])
                        weight[i] = 1 - tau;
                    else
                        weight[i] = tau;
            }
            
            iter = 1;
            VectorXd wold = weight;
            wold[1]++;
            while((wold - weight).norm() != 0 && iter < 50) {
                checkUserInterrupt();
                
                llt.compute(X.adjoint() * weight.asDiagonal() * X);
                betahat = llt.solve(X.adjoint() * weight.asDiagonal() * y);
                
                iter++;
                
                wold = weight;
                fitted = X*betahat;
                
                for(int i=0;i<n;i++)
                    if(y[i] < fitted[i])
                        weight[i] = 1 - tau;
                    else
                        weight[i] = tau;
                    
            }
            // # End Fit LAWS for given lambda
            
            double sig2 = (((weight.head(n).array())*((pow((((y.head(n) - fitted.head(n)).array())),2)).array())).sum() ) / (m - ((X*llt.solve(X.adjoint() * weight.asDiagonal())).diagonal().head(n)).sum()) ;
            
            VectorXd l0 = lambdashort;
            
            for(int k=0; k < glatterms.size(); k++) {
                checkUserInterrupt();
                
                MatrixXd partB;
                MatrixXd partDD;
                VectorXd partaa;
                
                if(center) {
                    // glatterms(k) -1; which glatterm to use; -1 because C++ starts vectors at 0
                    // NBPC() - 1; because NPBC is the starting point ion R; in C++ -1, because of starting vectors at 0
                    partB  = B.block(0,NBPC(glatterms(k)-1)-1,  B.rows(),NBP(glatterms(k)-1));
                    partDD = DDorig.block(NBPC(glatterms(k)-1)-1,NBPC(glatterms(k)-1)-1,NBP(glatterms(k)-1),NBP(glatterms(k)-1));
                    partaa = betahat.segment(NBPC(glatterms(k)-1)-1,NBP(glatterms(k)-1));
                    /*if(it == 1) {
                    Rcout << NBPC(glatterms(k)-1) << endl;
                    Rcout << NBP(glatterms(k)-1) << endl;
                    Rcout << partB.block(0,0,1,partB.cols()) << std::endl;
                    //Rcout << partDD << std::endl;
                    //Rcout << partaa << std::endl << std::endl << std::endl;
                    }*/
                }
                
                if(!center) {
                    MatrixXd B_red = B.block(0,1,B.rows(),B.cols());
                    MatrixXd DDorig_red = DDorig.block(1,1,DDorig.rows(),DDorig.cols());
                    VectorXd a_red = betahat.segment(1,n);
                    
                    partB  = B_red.block(0,NBPC(glatterms(k)-1),  B_red.rows(),NBP(glatterms(k)-1));
                    partDD = DDorig_red.block(NBPC(glatterms(k)-1)-1,NBPC(glatterms(k)-1)-1,NBP(glatterms(k)-1),NBP(glatterms(k)-1));
                    partaa = (a_red).segment(NBPC(glatterms(k)-1)-1,NBP(glatterms(k)-1));
                }
                
                if(NBP(glatterms(k)-1) != 0) {
                    VectorXd v = partDD * partaa;
                    
                    MatrixXd partX (partB.rows()+partDD.rows(),partB.cols());
                    
                    for(int i=0; i<partDD.cols();i++)
                        partDD.col(i) *= sqrt(lambdashort(glatterms(k)-1));
                    
                    partX << partB,partDD;
                    LLT<MatrixXd> partllt(AtA(partX));
                    
                    partllt.compute(partX.adjoint() * weight.asDiagonal() * partX);
                    double SumH = ((partX*partllt.solve(partX.adjoint() * weight.asDiagonal())).diagonal().head(n)).sum();
                    
                    double tau2 = (((pow(v.array(),2)).sum()) / SumH ) + 1e-6;
                    lambdashort(glatterms(k)-1) = std::max((sig2/tau2),1e-10);
                }
            }
            dc = (((((l0.array() + 1e-6).array()).log()/log(10) - ((lambdashort.array() + 1e-6).array()).log()/log(10)).abs()).array()).maxCoeff();
            dw = (weight-w0).norm();
            it++;
        }
        
        return List::create(Named("a") =betahat,
                            Named("lambdashort") = lambdashort,
                            Named("diag.hat.ma") = (X*llt.solve(X.adjoint() * weight.asDiagonal())).diagonal().head(n),
                            Named("fitted") = fitted.head(n),
                            Named("weight") = weight.head(n),
                            Named("it") = it,
                            Named("iter") = iter
        );
}

