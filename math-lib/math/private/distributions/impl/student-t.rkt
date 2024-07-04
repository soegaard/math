#lang typed/racket

(require "../../functions/incomplete-beta.rkt"
         "../../functions/beta.rkt"
         "../normal-dist.rkt"
         "../gamma-dist.rkt"         
         "../../../flonum.rkt"
         "normal-random.rkt")

(provide (rename-out [make-pdf         make-student-t-pdf]
                     [make-cdf         make-student-t-cdf]
                     [make-inverse-cdf make-student-t-inverse-cdf])
         sample-student-t)
                     

#;(require racket/performance-hint
         
         )

;; (require (only-in math/special-functions  beta beta-inc)
;;          (only-in math/flonum             flbracketed-root fl)
;;          (only-in math/distributions      sample normal-dist gamma-dist)
;;          (only-in math/statistics         mean variance)) ; for the tests

;;; Student t distribution

;; Parameters
;   μ   - location parameter
;   σ   - scale parameter
;   ν   - degrees of freedom

;; Domains
;   μ   - any real number
;   σ   - any positive real number
;   ν   - any positive real number


;; Probability Density Function (PDF)


; The density is proportional to
;             ν                           
; f(x) = ( -------- )^( (1+ν)/2 )         
;           (ν+x²)                        

; The proportionality constant is

;     Γ((ν+1)/2)                    
; -----------------                 
;  sqrt(πν) Γ(ν/2)                  

; Using the Β function this can also be written as

;          1                         
; --------------------               
;  sqrt(ν) B(1/2, ν/2)               

;; Reduction to one parameter

; If X ~ t(μ,σ,ν) then (X-μ)/σ ~ t(ν).
; This means we can concentrate on the one parameter version t(ν).
;     t(ν) = t(0, 1, ν)

;  t(ν)      - "the" Student t distribution
;  t(μ,σ,ν)  - the generalized Student t distribution
;              also called the location-scale-t-distribution


;; The Cumulative distribution function (CDF)

; For t>0
;   F(t) = 1 - 0.5 I_x(t) (ν/2, 0.5),
;   where x(t) = ν/(t²+ν).

; Here I is the regularized incomplete beta function.


;; Implementation Notes

; Names:  (student-t-dist              ν)       
;         (generalized-student-t-dist  μ σ ν)
; Note: (student-t-dist [μ 0] [σ 1] ν) could be implemented with case-lambda


; Type returned:  (ordered-dist Real Flonum)


; (beta-regularized z a b)
;    The regularized incomplete beta function I_z(a,b)
(: beta-regularized (-> Real Real Real Real))
(define (beta-regularized z a b) ; a, b > 0  
  ; The arguments for beta-inc
  ;   #f means integrate from 0
  ;   #t means regularized 
  (beta-inc a b z #f #t))
  
(: make-pdf : (case-> (Real           -> (-> Real Flonum))
                      (Real Real Real -> (-> Real Flonum))))
; TODO: check ν>0
(define make-pdf
  (case-lambda
    ; X ~ t(ν)
    [(ν)     (define proportionality-constant (/ 1. (* (sqrt ν) (beta 0.5 (/ ν 2.)))))
             (: pdf (-> Real Flonum))
             (define (pdf x)
               (define base (/ ν (+ ν (* x x))))
               (define expo (/ (+ 1. ν) 2.))         
               (fl (cast (* proportionality-constant (expt base expo)) Real)))
             pdf]
    ; Y ~ σX+μ
    [(μ σ ν) (define f (make-pdf ν))
             (λ (y)
               (define x (/ (- y μ) σ))
               (/ (f x) σ))]))


(: make-cdf : (case-> (Real           -> (-> Real Flonum))
                      (Real Real Real -> (-> Real Flonum))))
(define make-cdf
  (case-lambda
    ; X ~ t(ν)
    [(ν)     (define ν/2 (/ ν 2))
             (: sub (Real -> Real))
             (define (sub x) (/ ν (+ (* x x) ν)))
             (: cdf (Real -> Flonum))
             (define (cdf x)
               (cond
                 ; the distribution is symmetrical around x=0
                 [(= x 0.) 0.5]
                 ; reduce to the case x>0
                 [(< x 0.) (- 1. (cdf (- x)))]
                 ; general case
                 [else     (- 1. (* 0.5 (beta-regularized (sub x) ν/2 0.5)))]))
             cdf]
    ; Y ~ σX+μ
    [(μ σ ν) (define F (make-cdf ν))
             (λ (y)
               (define x (/ (- y μ) σ))
               (F x))]))


(: find-bracket : ((Real -> Real) Real Real  -> (values Real Real)))
(define (find-bracket h [a -1.] [b 1.] )
  ; Since the function h is monotone, this strategy works.
  (define ha (h a))
  (define hb (h b))
  (if (or (and (positive? ha) (negative? hb))
          (and (negative? ha) (positive? hb)))
      (values a b)
      (find-bracket h (* 2. a) (* 2. b))))

(: make-inverse-cdf : (case-> (Real           -> (-> Real Flonum))
                              (Real Real Real -> (-> Real Flonum))))
(define make-inverse-cdf
  (case-lambda
    ; X ~ t(ν)
    [(ν)     (case ν
               ; special cases
               [(1 2 4)   (: inv-F : (Real -> Flonum))
                          (define inv-F
                           (case ν
                             [(1) (λ ([p : Real])
                                    (fl (cast (tan (* pi (- p 0.5))) Real)))]
                             [(2) (λ ([p : Real])
                                    (define α (* 4. p (- 1. p)))
                                    (fl (cast (* 2. (- p 0.5) (sqrt (/ 2 α))) Real)))]
                             [(4) (λ ([p : Real])
                                    (define α (* 4. p (- 1. p)))
                                    (define q (/ (cos (/ (acos (sqrt α)) 3.)) (sqrt α)))
                                    (fl (cast (* (sgn (- p 0.5)) 2. (sqrt (- q 1.))) Real)))]
                             [else (λ ([p : Real]) 0.0)])) ; happy type checking
                         (λ ([p : Real])
                           (case p
                             [(0) -inf.0]
                             [(1) +inf.0]
                             [else (inv-F p)]))]
               ; general
               [else (define F (make-cdf ν))
                     (λ (p)
                       (cond
                         [(= p 0) -inf.0]
                         [(= p 1) +inf.0]
                         [else
                          ; a root of F(x)=p is the inverse of F in p
                          (: h : (Real -> Flonum))
                          (define (h x) (fl (- (F x) p)))
                          ; find interval in which the root lies
                          (define-values (a b) (find-bracket h -1. 1.))
                          ; find the root
                          (flbracketed-root h (fl a) (fl b))]))])]
    ; Y ~ σX+μ
    [(μ σ ν) (define inv-F (make-inverse-cdf ν))
             (λ (p)
               (define x (inv-F p))
               (define y (+ (* σ x) μ))
               (fl y))]))

(: sample-student-t : (case-> (Real           Integer -> FlVector)
                              (Real Real Real Integer -> FlVector)))
(define sample-student-t
  (case-lambda
    ; X ~ t(ν)
    
    [(ν n) (cond
             [(n . < . 0)  (raise-argument-error 'sample-student-t "Natural" 1 n)]
             [else    

              (define ν/2 (fl (/ ν 2.)))
              ; Note: Our gamma distribution has a shape parameter.
              ;       A shape parameter of 2 corresponds to a a rate of 1/2.
              (define Xs  (flnormal-sample 0. 1. n))
              (define X²s (flgamma-sample ν/2 2. n))
              (build-flvector n
                              (λ (i)
                                (define X  (flvector-ref Xs  i))
                                (define X² (flvector-ref X²s i))
                                (fl (cast (/ X (sqrt (/ X² ν))) Real))))])]

    ; Y ~ σX+μ
    [(μ σ ν n) (cond
                 [(n . < . 0)  (raise-argument-error 'sample-student-t "Natural" 3 n)]
                 [else    
                  (define ν/2 (fl (/ ν 2.)))
                  (define Xs  (flnormal-sample 0. 1. n))
                  (define X²s (flgamma-sample ν/2 2. n))
                  (build-flvector n
                                  (λ (i)
                                    (define X  (flvector-ref Xs  i))
                                    (define X² (flvector-ref X²s i))
                                    (define x (/ X (sqrt (/ X² ν))))
                                    (fl (cast (+ (* σ x) μ) Real))))])]))

;;;
;;; Tests
;;;


(: nearly-equal? : (Real Real Real -> Boolean))
(define (nearly-equal? eps x y)
  (<= (abs (- x y)) eps))

;; Numerical test cases were computed by the free `wolframscript`.

;; All tests are expected to return $t.
;; If a set of tests results in #f, change `and` to `list` to see the
;; individual results.

#;(list "Density - PDF"
      (and
       ; N[PDF[StudentTDistribution[1], 0], 30]
       (nearly-equal? (expt 2 -54) ((make-pdf 1) 0)      0.3183098861837907)
       (nearly-equal? (expt 2 -55) ((make-pdf 2) 1)      (/ 1 (* 3 (sqrt 3))))
       ; generalized
       (nearly-equal? (expt 2 -55) ((make-pdf 2 2 1) 0)  0.07957747154594767)
       (nearly-equal? (expt 2 -55) ((make-pdf 2 2 1) 1)  0.12732395447351627)
       (nearly-equal? (expt 2 -55) ((make-pdf 3 4 5) 0)  0.06892452901798418)
       (nearly-equal? (expt 2 -55) ((make-pdf 3 4 5) 1)  0.08197963283068663))
      "Cumulative - CDF"
      (and
       ; N[CDF[StudentTDistribution[1], 0], 30]
       (equal?                      ((make-cdf 1)  0)      0.5)
       (equal?                      ((make-cdf 1)  1)      0.75)
       (equal?                      ((make-cdf 1) -1)      0.25)
       (equal?                      ((make-cdf 2)  0)      0.5)
       (equal?                      ((make-cdf 2)  1)      0.7886751345948129)
       (nearly-equal? (expt 2 -55)  ((make-cdf 2) -1)      0.2113248654051871)
       ; generalized
       (nearly-equal? (expt 2 -55)  ((make-cdf 1 2 3) -1)  0.19550110947788532)
       (nearly-equal? (expt 2 -55)  ((make-cdf 1 2 3)  0)  0.3257239824240755)
       (nearly-equal? (expt 2 -55)  ((make-cdf 1 2 3)  1)  0.5))
      " Inverse Cumulative - Inverse CDF"
      ; Example to get expected result:
      ;   N[InverseCDF[StudentTDistribution[2], 1/10], 30]
      (and
        ; Special case ν=1
        (equal?                           ((make-inverse-cdf 1)  0)   -inf.0)
        (equal?                           ((make-inverse-cdf 1)  1)   +inf.0)
        (equal?                           ((make-inverse-cdf 1)  0.5) 0.) 
        (nearly-equal? (expt 2 -51)       ((make-inverse-cdf 1)  0.1) -3.0776835371752536)
        (nearly-equal? (expt 2 -51)       ((make-inverse-cdf 1)  0.9)  3.0776835371752536)
        ; Special case ν=2
        (equal?                           ((make-inverse-cdf 2)  0)   -inf.0)
        (equal?                           ((make-inverse-cdf 2)  1)   +inf.0)
        (equal?                           ((make-inverse-cdf 2)  0.5) 0.) 
        (nearly-equal? (expt 2 -51)       ((make-inverse-cdf 2)  0.1) -1.8856180831641267)
        (nearly-equal? (expt 2 -51)       ((make-inverse-cdf 2)  0.9)  1.8856180831641267)
        ; Special case ν=4
        (equal?                           ((make-inverse-cdf 4)  0)   -inf.0)
        (equal?                           ((make-inverse-cdf 4)  1)   +inf.0)
        (equal?                           ((make-inverse-cdf 4)  0.5) 0.) 
        (nearly-equal? (expt 2 -51)       ((make-inverse-cdf 4)  0.1) -1.5332062740589438)
        (nearly-equal? (expt 2 -51)       ((make-inverse-cdf 4)  0.9)  1.5332062740589438)
        ; General case
        (equal?                           ((make-inverse-cdf 3)  0.5) 0.) 
        (nearly-equal? (expt 2 -52)       ((make-inverse-cdf 3)  0.1) -1.6377443536962102)
        (nearly-equal? (expt 2 -52)       ((make-inverse-cdf 3)  0.9)  1.6377443536962102)
        (equal?                           ((make-inverse-cdf 5)  0)   -inf.0)
        (equal?                           ((make-inverse-cdf 5)  1)   +inf.0)
        (equal?                           ((make-inverse-cdf 5)  0.5) 0.) 
        (nearly-equal? (expt 2 -52)       ((make-inverse-cdf 5)  0.1) -1.475884048824481)
        (nearly-equal? (expt 2 -52)       ((make-inverse-cdf 5)  0.9)  1.475884048824481)
        ; Three parameters
        (nearly-equal? (expt 2 -55)  ((make-inverse-cdf 1 2 3)  0.5)  1)
        (nearly-equal? (expt 2 -51)  ((make-inverse-cdf 1 2 3)  0.1)  -2.2754887073924204)
        (nearly-equal? (expt 2 -50)  ((make-inverse-cdf 1 2 3)  0.9)   4.27548870739242)))
