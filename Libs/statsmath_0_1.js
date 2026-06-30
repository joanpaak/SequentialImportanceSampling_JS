/**
 * Collection of r-style functions for common probability
 * density and mass functions.
 */

class Stats{
    static dunif(x, min, max){
        if(x < min | x > max) return 0;

        return 1.0 / (max - min);
    }

    static runif(n, min, max){
        let x = new Array(n);

        for(let i = 0; i < n; i++){
            x[i] = Math.random() * (max - min) + min;
        }

        return x;
    }

    //static dlogis(x){
    // TODO: implement
    //}

    // https://en.wikipedia.org/wiki/Logistic_distribution
    static plogis(p, mu, s){
        return 1.0 / (1 + Math.exp(-(p - mu) / s));
    }

    // https://en.wikipedia.org/wiki/Logistic_distribution
    static qlogis(q, mu, s){
        return mu + s * Math.log(q / (1 - q));
    }

    // https://en.wikipedia.org/wiki/Beta_distribution
    static dbeta(x, alpha, beta){
        if(x < 0 | x > 1) return 0;

        let normConst = Functions.beta(alpha, beta);

        return Math.pow(x, alpha - 1) * Math.pow(1 - x, beta - 1) / normConst;
    }

    // https://en.wikipedia.org/wiki/Beta_distribution#Random_variate_generation
    static rbeta(n, alpha, beta){
        let x = Stats.rgamma(n, alpha, 1);
        let y = Stats.rgamma(n, beta, 1);
        // 
        let z = new Array(n);

        for(let i = 0; i < n; i++){
            z[i] = x[i] / (x[i] + y[i]);
        }

        return z;
    }

    // qbeta TODO
    // pbeta TODO

    // https://en.wikipedia.org/wiki/Normal_distribution
    static dnorm(x, mu, sigma){
        let nc = 1.0 / Math.sqrt(2 * Math.PI * (sigma * sigma));
        let p = Math.exp(-((x - mu) * (x - mu)) / (2 * (sigma * sigma)));

        return nc * p;
    }

    static qnorm(p, mu, sigma){
        return Stats.invPhi(p) * sigma + mu;
    }

    // The Shore (1982) algorith taken from: https://en.wikipedia.org/wiki/Normal_distribution
    static invPhi(p){
        if(p >= 0.5){
            return 5.5556 * (1.0 - Math.pow((1.0 - p) / p, 0.1186));
        }

        return(-Stats.invPhi(1 - p));
    }

    static pnorm(x, mu, sigma){
        return Stats.phi((x - mu) / sigma);
    }

    // Algorithm by Zelen & Severo (1964), from 
    // https://en.wikipedia.org/wiki/Normal_distribution 
    static phi(x){
        let t, p;

        const b0 = 0.2316419;
        const b1 = 0.319381530;
        const b2 = -0.356563782;
        const b3 = 1.781477937;
        const b4 = -1.821255978;
        const b5 = 1.330274429;

        if(x >= 0){
            t = 1.0 / (1 + b0 * x);
            p = 1 - Stats.dnorm(x, 0, 1) * 
                (b1 * t + b2 * (t * t) + 
                b3 * (t * t * t) + 
                b4 * (t * t * t * t) + 
                b5 * (t * t * t * t * t));
            
            return p;
        }
        
        return 1 - Stats.phi(-x);
    }

    static rnorm(n, mu, sigma){
        let x = new Array(n);

        for(let i = 0; i < n; i++){
            x[i] = Stats.invPhi(Math.random()) * sigma + mu;
        }

        return x;        
    }

    // https://en.wikipedia.org/wiki/Log-normal_distribution
    static dlnorm(x, mu, sigma){
        let normConst = 1.0 / (x * sigma * Math.sqrt(2 * Math.PI));
        let p = Math.exp(- Math.pow(Math.log(x) - mu, 2) / (2 * sigma * sigma));

        return normConst * p;
    }

    static rlnorm(n, mu, sigma){
        // TODO: Just use the normal generator and exponentiate samples
    }

    // https://en.wikipedia.org/wiki/Student%27s_t-distribution
    static dt(x, mu, s, nu){
        x = (x - mu) / s;

        let normConst = Functions.gamma((nu + 1) / 2) / (Math.sqrt(Math.PI * nu) * 
            Functions.gamma(nu / 2));
        let p = Math.pow(1 + ((x * x) / nu), -((nu + 1)/2));

        return normConst * p;
    }

    // NOTE: following R stan, scale parameter s does not correspond
    // to standard deviation; expect for most cases standard deviation 
    // to be slightly more than s.
    // From:
    // Polar generation of random variates withe the t-distribution
    // Ralph W. Bailey (1994). Mathematics of Compuation, 62(206), 779 - 781
    static rt(n, mu, s, nu){
        let U, V, C, R, sigma, x;

        x = new Array(n);

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
     * https://en.wikipedia.org/wiki/Gamma_distribution
     * 
     * NOTE: R uses shape and rate parameterisation, in which mean
     * of the deviates is shape / rate
     * @param {Number} x 
     * @param {*} shape 
     * @param {*} rate 
     * @returns 
     */
    static dgamma(x, shape, rate){
        if(x < 0) return 0;

        let normConst = Math.pow(rate, shape) / Functions.gamma(shape);
        let p = Math.pow(x, shape - 1) * Math.exp(- rate * x);

        return normConst * p;
    }

    /**
     * From: Marsaglia & Tsang (2000). A simple method for generating 
     *   gamma variables.
     * @param {Number} shape 
     * @returns x ~ gamma(shape, 1);
     */
    static gammaDeviate(shape){
        let x, c, d, g, i, s, v;
        i = 0;
        d = shape - (1/3);
        c = 1.0 / Math.sqrt(9 * d);
        while(true){
            x = this.rnorm(1, 0, 1);
            v = (1 + c * x);
            v = v * v * v;

            if(v > 0){
                g = 0.5 * (x * x) + d - d * v + d * Math.log(v);
                s = Math.log(Math.random());

                if(s < g){
                    return d * v
                }
            }

            i++;
            if(i > 10000){
                throw new Error("Gamma generator is stuck!");
            }
        }
    }

    /**
     * 
     */
    static rgamma(n, shape, rate){
        let x = new Array(n);

        for(let i = 0; i < n; i++){
            x[i] = Stats.gammaDeviate(shape) / rate;
        }

        return x;
    }

    // qgamma TOOD
    // pgamma TODO

    // https://en.wikipedia.org/wiki/Dirichlet_distribution
    static ddirichlet(x, alpha){
        function B(alpha){
            let prod = 1.0;
            let sum = 0;

            for(let i = 0; i < alpha.length; i++){
                prod *= Functions.gamma(alpha[i]);
                sum += alpha[i];
            }

            return prod / Functions.gamma(sum);
        }

        let normConst = 1.0 / B(alpha);
        let p = 1.0;

        for(let i = 0; i < alpha.length; i++){
            p *= Math.pow(x[i], alpha[i] - 1);
        }

        return normConst * p;
    }

    // https://en.wikipedia.org/wiki/Dirichlet_distribution#Random_variate_generation
    static rdirichlet(n, alpha){
        let x = new Array(n);

        for(let j = 0; j < n; j++){
            let y = new Array(alpha.length);
            let sum = 0;
    
            for(let i = 0; i < alpha.length; i++){
                y[i] = Stats.rgamma(1, alpha[i], 1)[0];
                sum += y[i];
            }
    
            for(let i = 0; i < alpha.length; i++){
                y[i] /= sum;
            }

            x[j] = y;
        }

        return x;
    }

    /**
     * 
     * @param {Number} x 
     * @param {Array} alpha 
     * @returns 
     */
    static dcategorical(x, alpha){
        return alpha[x];
    }

    static rcategorical(n, alpha){
        let x = new Array(n);
        let q = Functions.cumsum(alpha);
        
        for(let i = 0;i < n; i++){
            let s = Math.random();
            let ind = 0;
    
            while(q[ind] < s){
                ind++;
                if(ind > alpha.length){
                    throw new Error("Error in rcategorical: varible ind > alpha.length");
                }
            } 

            x[i] = ind;
        }

        return x;
    }

    static rbinom(numSamples, n,  p){
        let d = new Array(n + 1);
        let x = new Array(numSamples);

        for(let i = 0; i < (n +1); i++){
            d[i] = Stats.dbinom(i, n, p);
        }

        x = Stats.rcategorical(numSamples, d);

        return x;
    }

    // https://en.wikipedia.org/wiki/Binomial_distribution
    static dbinom(x, n, p){ 
        return Functions.choose(n, x) * 
            (Math.pow(p, x) * Math.pow(1.0 - p, n - x));
    }

    static rbernoulli(n, p){
        let x = new Array(n);
        
        for(let i = 0; i < n; i++){
            if(Math.random() <= p){
                x[i] = 1;
            } else {
                x[i] = 0;
            }
        }

        return x;
    }

    static dbernoulli(x, p){
        return x * p + (1 - x) * (1 - p);
    }
}

class Functions{
    static cumsum(x){
        let q = new Array(x.length);
        q[0] = x[0];

        for(let i = 1; i < x.length; i++){
            q[i] = q[i - 1] + x[i];
        }

        return q;
    }


    // https://en.wikipedia.org/wiki/Beta_distribution
    static beta(a, b){
        return (Functions.gamma(a) * Functions.gamma(b)) / Functions.gamma(a + b);
    }

    // https://en.wikipedia.org/wiki/Stirling%27s_approximation#Versions_suitable_for_calculators
    static gamma(x){
        return Math.sqrt((2 * Math.PI) / x) * 
            Math.pow((1.0 / Math.E) * 
            (x + 1.0 / (12 * x - 1.0 / (10 * x))) , x); 
    }

    static mean(x){
        let sum = 0;

        for(let i = 0; i < x.length; i++){
            sum += x[i];
        }

        return sum / x.length;
    }
    
    static var(x){
        let mu = this.mean(x);
        let sum = 0;

        for(let i = 0; i < x.length; i++){
            sum += (x[i] - mu) * (x[i] - mu); 
        }

        return sum / x.length;
    }

    static sd(x){
        return Math.sqrt(this.var(x));
    }


    // n! / (k! * (n -k)!)
    static choose(n, k){
        return this.factorial(n) / (this.factorial(k) * this.factorial(n - k));
    }

    static factorial(x){
        if(x === 0) return 1;
        if(x === 1) return 1;
        if(x === 2) return 2;
        let p = 2;

        for(let i = 3; i <= x; i++){
            p *= i;
        }

        return p;
    }
}
