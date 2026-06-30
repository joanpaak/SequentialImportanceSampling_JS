/**
 * SEQUENTIAL IMPORTANCE SAMPLING WITH TEMPERING
 * version 0.1
 */

class Matrix{
    /**
     * You can initialize this class in two ways:
     *   1) give it an array of arrays which is interpreted as matrix
     *   2) give it an object with properties ncol and nrow which are
     *      used for initializing an empty matrix.
     * @param {} x 
     */
    constructor(x){
        if(Array.isArray(x)){
            this.x = x;
            this.nrow = this.x.length;
            this.ncol = this.x[0].length;
        } else {
            this.ncol = x.ncol;
            this.nrow = x.nrow;

            this.x = new Array(this.nrow);

            for(let i = 0; i < this.nrow; i++){
                this.x[i] = new Array(this.ncol);
            }
        }
    }

    get(row, col){
        return this.x[row][col];
    }

    set(row, col, x){
        this.x[row][col] = x;
    }

    getData(){
        return this.x;
    }

    getCol(ind){
        let c = new Array(this.nrow);

        for(let i = 0; i < this.nrow; i++){
            c[i] = this.x[i][ind];
        }

        return c;
    }

    getRow(ind){
        return this.x[ind];
    }

    setCol(ind, x){
        if(!Array.isArray(x)) x = [x];

        for(let i = 0; i < this.nrow; i++){
            this.x[i][ind] = x[i];
        }
    }

    setRow(ind, x){
        if(!Array.isArray(x)) x = [x];

        for(let i = 0; i < x.length; i++){
            this.x[ind][i] = x[i];
        }
    }

    appendRow(row){
        if(row.length != this.ncol){
            throw new Error("Cannot append row: dimension mismatch");
        }

        this.nrow++;
        this.x[this.nrow - 1] = row;
    }

    appendCol(col){
        if(col.length != this.nrow){
            throw new Error("Cannot append col: dimension mismatch");
        } 

        for(let i = 0; i < this.nrow; i++){
            this.x[i][this.ncol] = col[i];
        }
        
        this.ncol++;
    }
}

class SIS{
    /**
     * 
     * @param {Object} settings
     * @param {Function} settings.likelihood
     * @param {Function} settings.prior
     * @param {Function} settings.drawFromPrior
     * @param {Number} settings.nParticles 
     * @param {Object} opt 
     * @param {Number} opt.nu 
     * @param {Number} opt.rejuvenationLimit
     * @param {Number} opt.temperingLimit
     * @param {Number} opt.k
     * @param {Boolean} opt.autoAdjustK
     * @param {Number} opt.minK
     * @param {Number} opt.maxK
     * @param {Number} opt.logging
     */
    constructor(settings, opt){
        if(settings.nParticles === undefined){
            throw new Error("You must specify sample size");
        }

        // TODO: Mayhaps tests for seeing if these are functions
        if(settings.drawFromPrior === undefined){
            throw new Error("You must specify function for sampling from prior");
        }

        if(settings.prior === undefined){
            throw new Error("You must specify function for calculating prior probabilities");
        }

        if(settings.likelihood === undefined){
            throw new Error("You must specify function for calculating likelihoods");
        }

        //

        this.drawFromPrior = settings.drawFromPrior;
        this.prior = settings.prior;
        this.likelihood = settings.likelihood;

        /**
         * 
         * @param {Matrix} theta 
         * @param {Matrix} y 
         * @returns Array of posterior probabilities
         */
        this.posterior = (theta, y) => {
            let postProb = this.prior(theta);

            for(let i = 0; i < y.nrow; i++){
                postProb = this.arrMult(postProb, 
                    this.likelihood(y.getRow(i), theta));
            }

            this.checkAndFixProbArray(postProb);

            return postProb;
        }

        this.nParticles = settings.nParticles;
        this.y = null;
        // NOTE: 1st dim = rows, 2nd = cols
        this.theta = this.drawFromPrior(this.nParticles);
        this.nDim = this.theta.ncol;
        this.w = new Array(this.nParticles);

        for(let i = 0; i < this.nParticles; i++){
            this.w[i] = 1.0 / this.nParticles;
        }

        this.opt = {
            nu : 10,
            temperingLimit : 0.50,
            rejuvenationLimit : 0.75,
            k : 5,
            minK : 2,
            maxK : 40,
            nRejuvenationSteps : 1,
            autoAdjustK : true,
            logging : false
        }
   
        if(opt != undefined){
            Object.keys(opt).forEach(key => {
                this.opt[key] = opt[key];
            });
        }

        this.temperingLimitN = this.nParticles * this.opt.temperingLimit;
        this.rejuvenationLimitN = this.nParticles * this.opt.rejuvenationLimit;
        
        this.log = {
            k : [["N", "k"]],
            nEff : [["N", "N eff"]],
            nAccepted : [["N", "N Accepted"]],
            particleSets : new Array()
        }

        this.warnings = []; // TODO
    }

    /**
     * 
     * @param {Array} y 
     * @param {Object} opt
     * @param {Boolean} opt.forceRejuvenation
     * @param {Boolean} opt.forceTempering 
     */
    addObservation(y, opt){
        // TODO THIS IS HORRIBLE!!!!! HORRIBLE!!!!
        if(this.y === null){
            if(Array.isArray(y)){
                this.y = new Matrix({nrow : 1, ncol : y.length});
                this.y.setRow(0, y);
            } else {
                this.y = y;
            }
        } else {
            if(Array.isArray(y)){
                this.y.appendRow(y);
            } else {
                this.y.appendRow(y.getRow(0));
            }
        }

        let p = this.likelihood(y, this.theta);
        let wUpdated = this.arrMult(p, this.w);
        wUpdated = this.normArr(wUpdated);

        let nEff = this.nEff(wUpdated);

        let nEffAboveRejuvenation = nEff > this.rejuvenationLimitN;
        let nEffAboveTempering = nEff > this.temperingLimitN;

        if(nEffAboveRejuvenation & nEffAboveTempering){
            this.w = wUpdated;
        } else if(nEffAboveRejuvenation & !nEffAboveTempering){
            this.doTempering(y);
        } else if(!nEffAboveRejuvenation & nEffAboveTempering){
            for(let j = 0; j < this.opt.nRejuvenationSteps; j++){
                this.rejuvenate();
            }
        } else if(!nEffAboveRejuvenation & !nEffAboveTempering){
            this.doTempering(y, true);
        }

        this.checkAndFixProbArray(this.w);

        if(this.opt.logging){
            this.log.nEff.push(
                [
                    this.y.nrow,
                    this.nEff(this.w)
                ]
            );
        }
    }

    adjustK(y){
        let lh = this.likelihood(y, this.theta);
        let temperedLh = new Array(this.nParticles);

        for(let k = this.opt.minK; k <= this.opt.maxK; k++){
            for(let i = 0; i < this.nParticles; i++){
                temperedLh[i] = Math.pow(lh[i], 1.0 / k);
            }

            let wUpdated = this.arrMult(temperedLh, this.w);
            wUpdated = this.normArr(wUpdated);

            let nEff = this.nEff(wUpdated);

            if(nEff > this.opt.temperingLimit){
                this.opt.k = k;
                return;
            }
        }

        this.warnings.push("Reached maxK when adjusting K");
        this.opt.k = this.opt.maxK;

        return;
    }

    /**
     * 
     * @param {Array} y 
     * @param {Boolean} forceRejuvenation 
     */
    doTempering(y, forceRejuvenation){
        let lh;

        if(this.opt.autoAdjustK){
            this.adjustK(y);
        }

        if(this.opt.logging){
            this.log.k.push(
                [this.y.nrow, this.opt.k]
            );
        }

        for(let i = 0; i < this.opt.k; i++){
            lh = this.likelihood(y, this.theta);
            
            for(let j = 0; j < this.nParticles; j++){
                lh[j] = Math.pow(lh[j], 1.0 / this.opt.k);
            }

            this.w = this.arrMult(this.w, lh);
            this.w = this.normArr(this.w);
            this.checkAndFixProbArray(this.w);

            if(this.nEff(this.w) < this.rejuvenationLimitN | forceRejuvenation){
                for(let j = 0; j < this.opt.nRejuvenationSteps; j++){
                    this.rejuvenate();
                }
            }
        }
    }

    rejuvenate(){
        let pprop;
        let mhrat;
        let nAccepted = 0;
        let pcur;
        let margMus = new Array(this.nDim);
        let margSDs = new Array(this.nDim);

        for(let i = 0; i < this.nDim; i++){
            margMus[i] = this.mean(this.theta.getCol(i), this.w);
            margSDs[i] = this.sd(this.theta.getCol(i), this.w);
        }

        let rsIndices = this.sampleMultinomial(this.nParticles, this.w);

        let thetaRs = new Matrix({
            ncol : this.nDim,
            nrow : this.nParticles
        });
        
        for(let i = 0; i < rsIndices.length; i++){
            thetaRs.setRow(i, this.theta.getRow(rsIndices[i]));
        }

        this.theta = thetaRs;

        let proposals = new Matrix({
            ncol : this.nDim,
            nrow : this.nParticles
        });

        for(let i = 0; i < this.nDim; i++){
            proposals.setCol(i, this.tRandom(
                this.nParticles, margMus[i], margSDs[i], this.opt.nu));
        }

        let pCurrent = this.posterior(this.theta, this.y);
        let pProposal = this.posterior(proposals, this.y);

        for(let i = 0; i < this.nParticles; i++){
            mhrat = this.multiTPDF(this.theta.getRow(i), margMus, margSDs, this.opt.nu) /
                this.multiTPDF(proposals.getRow(i), margMus, margSDs, this.opt.nu);

            let pAccept = (pProposal[i] / pCurrent[i]) * mhrat;
            if(isNaN(pAccept)) pAccept = 0;
                
            // TODO: Is it more efficient to do this in its own loop?
            if(Math.random() < pAccept){
                this.theta.setRow(i, proposals.getRow(i));
                nAccepted++;
            }
        }

        for(let i = 0; i < this.nParticles; i++){
            this.w[i] = 1.0 / this.nParticles;
        }

        if(this.opt.logging){
            this.log.nAccepted.push(
                [this.y.nrow, nAccepted]
            );
        }
    }

    getMarginalMus(){
        let margMus = new Array(this.nDim);
        
        for(let i = 0; i < this.nDim; i++){
            margMus[i] = this.mean(this.theta.getCol(i), this.w);
         }

         return margMus;
    }

    getMarginalSDs(){
        let margSDs = new Array(this.nDim);

        for(let i = 0; i < this.nDim; i++){
            margSDs[i] = this.sd(this.theta.getCol(i), this.w);
        }

        return margSDs;
    }

    /**
     * Samples elements from the theta matrix according to their current weights.
     * @returns 
     */
    getIIDSample(){
        let inds = this.sampleMultinomial(this.nParticles, this.w);
        let theta = new Matrix({ncol : this.nDim, nrow : this.nParticles});

        for(let i = 0; i < this.nParticles; i++){
            theta.setRow(i, this.theta.getRow(inds[i]));
        }

        return theta;
    }

    /* AUXILIARY STUFF */

    /**
     * Mean of array x with weights p
     * @param {Array} x 
     * @param {Array} p 
     * @returns 
     */
    mean(x, p){
        let a = 0;

        for(let i = 0; i < x.length; i++){
            a += x[i] * p[i];
        }

        return a;
    }

    /**
     * Standard deviation of an array x with weights p
     * @param {Array} x 
     * @param {Array} p 
     * @returns 
     */
    sd(x, p){
        let m = this.mean(x, p);
        let s = 0;

        for(let i = 0; i < x.length; i++){
            s += ((x[i] - m) * (x[i] - m)) * p[i];
        }

        return Math.sqrt(s);
    }

    /**
     * Sample n deviates from distribution p
     * @param {Number} n 
     * @param {Array} p 
     * @returns 
     */
    sampleMultinomial(n, p){
        let q = new Array(p.length);
        let u = new Array(n);
        let s = new Array(n);

        q[0] = p[0];

        for(let i = 1; i < p.length; i++){
            q[i] = q[i - 1] + p[i];
        }

        for(let i = 0; i < n; i++){
            u[i] = Math.random();
        }

        for(let i = 0; i < n; i++){
            for(let j = 0; j < q.length; j++){
                if(u[i] < q[j]){
                    s[i] = j;
                    break;
                }
            }
        }

        return s;
    }

    /**
     * 
     * @param {Number} n 
     * @param {Number} mu 
     * @param {Number} sigma 
     * @param {Number} nu 
     * @returns 
     */
    tRandom(n, mu, s, nu){
        let U, V, C, R;

        let x = new Array(n);

        for(let i = 0; i < n; i++){
            let W = 2;
    
            while(W > 1){
                U = 2 * Math.random() - 1;
                V = 2 * Math.random() - 1;
                W = U * U + V * V;
            }
    
            C = U / Math.sqrt(W);
            R = Math.sqrt(nu * (Math.pow(W, -2/nu) - 1));

            x[i] = (R * C) * s + mu;
        }

        return x;
    }

    /**
     * Unnormalized density of the t-distribution
     * @param {Number} x 
     * @param {Number} mu 
     * @param {Number} s
     * @param {Number} nu 
     */
    tPDF(x, mu, s, nu){
        x = (x - mu) / s;

        let d = Math.pow(1 + (x * x) / nu, -(nu + 1) / 2);

        return d;
    }

    /**
     * Unnormalized PDF of the multivariate t-distribution, without
     * not taking correlations into account.
     * @param {Array} x 
     * @param {Array} mu 
     * @param {Array} s 
     * @param {Number} nu 
     * @returns 
     */
    multiTPDF(x, mu, s, nu){
        let d = 1.0;

        for(let i = 0; i < x.length; i++){
            d *= this.tPDF(x[i], mu[i], s[i], nu);
        }

        return d;
    }

    /**
     * Estimates effective sample size from an array
     * of normalized weights.
     * @param {Array} w 
     */
    nEff(w){
        let ss = 0;

        for(let i = 0; i < w.length; i++){
            ss += (w[i] * w[i]);
        }

        return 1.0 / ss;
    }

    /**
     * Normalize array so that sum(a) = 1.0
     * @param {Array} a 
     */
    normArr(a){
        let s = 0;
        let a2 = new Array(a.length);

        for(let i = 0; i < a.length; i++){
            s += a[i];
        }

        for(let i = 0; i < a.length; i++){
            a2[i] = a[i] / s;
        }

        return a2;
    }

    /**
     * Element-wise multiplication of arrays
     * @param {Array} a1 
     * @param {Array} a2 
     */
    arrMult(a1, a2){
        let a3 = new Array(a1.length);

        for(let i = 0; i < a1.length; i++){
            a3[i] = a1[i] * a2[i];
        }

        return a3;
    }

    /**
     * 
     * @param {Array} a 
     */
    checkAndFixProbArray(a){
        for(let i = 0; i < a.length; i++){
            if(isNaN(a[i]) | a[i] < 0){
                this.warnings.push(a[i] + " in array");
                a[i] = 0;
            }
        }
    }
}

