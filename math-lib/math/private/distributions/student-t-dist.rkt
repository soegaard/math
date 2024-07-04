#lang typed/racket/base

(require racket/performance-hint
         racket/promise
         "../../flonum.rkt"
         "../unsafe.rkt"
         "impl/normal-pdf.rkt"
         "impl/normal-cdf.rkt"
         "impl/normal-inv-cdf.rkt"
         "impl/normal-random.rkt"
         "impl/student-t.rkt"
         "dist-struct.rkt"
         "utils.rkt")

(provide Student-T-Dist
         student-t-dist
         ; student-t-dist-mean
         ; student-t-dist-stddev
         )




;; ===================================================================================================
;; Distribution object

(define-real-dist: student-t-dist Student-T-Dist
  student-t-dist-struct ([location : Flonum] [scale : Flonum] [freedom : Flonum] )) ; μ, σ, ν

(begin-encourage-inline
  (: student-t-dist (case-> (          -> Student-T-Dist)
                            (Real      -> Student-T-Dist)
                            (Real Real -> Student-T-Dist)))
  (define student-t-dist
    (case-lambda
      [(ν) (let ([ν (fl ν)])
             (define pdf     (make-student-t-pdf ν))
             (define cdf     (make-student-t-cdf ν))
             (define inv-cdf (make-student-t-inverse-cdf ν))
             (define sample  (case-lambda:
                               [()               (unsafe-flvector-ref (sample-student-t ν 1) 0)]
                               [([n : Integer])  (flvector->list (sample-student-t ν n))]))
             (student-t-dist-struct pdf sample cdf inv-cdf -inf.0 +inf.0 (delay ν) 0 1 ν))]
      [(μ σ ν) (let ([ν (fl ν)])
                 (define pdf     (make-student-t-pdf μ σ ν))
                 (define cdf     (make-student-t-cdf μ σ ν))
                 (define inv-cdf (make-student-t-inverse-cdf μ σ ν))
                 (define sample  (case-lambda:
                                   [()               (unsafe-flvector-ref (sample-student-t μ σ ν 1) 0)]
                                   [([n : Integer])  (flvector->list (sample-student-t μ σ ν n))]))
                 (student-t-dist-struct pdf sample cdf inv-cdf -inf.0 +inf.0 μ σ ν))])))
