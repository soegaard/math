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
         student-t-dist-location   ; μ
         student-t-dist-scale      ; σ
         student-t-dist-freedom)   ; ν


;; ===================================================================================================
;; Distribution object

(define-real-dist: student-t-dist Student-T-Dist
  student-t-dist-struct ([location : Real] [scale : Real] [freedom : Real] )) ; μ, σ, ν

(begin-encourage-inline
  (: student-t-dist (case-> (Real           -> Student-T-Dist)
                            (Real Real Real -> Student-T-Dist)))
  (define student-t-dist
    (case-lambda
      [(ν)     (let ([ν (fl ν)])
                 (cond
                   [(< ν 0) (raise-argument-error 'student-t-dist "Positive Real (ν)" 0 ν)]
                   [else    (define pdf     (make-student-t-pdf ν))
                            (define cdf     (make-student-t-cdf ν))
                            (define inv-cdf (make-student-t-inverse-cdf ν))
                            (define sample  (case-lambda:
                                              [()               (unsafe-flvector-ref (flstudent-t-sample ν 1) 0)]
                                              [([n : Integer])  (flvector->list (flstudent-t-sample ν n))]))
                            (define median (delay ν))
                            (student-t-dist-struct pdf sample cdf inv-cdf -inf.0 +inf.0 median 0 1 ν)]))]
      [(μ σ ν) (let ([μ (fl μ)] [σ (fl σ)] [ν (fl ν)])
                 (cond
                   [(< ν 0) (raise-argument-error 'student-t-dist "Positive Real (ν)" 2 μ σ ν)]
                   [else    (define pdf     (make-student-t-pdf μ σ ν))
                            (define cdf     (make-student-t-cdf μ σ ν))
                            (define inv-cdf (make-student-t-inverse-cdf μ σ ν))
                            (define sample  (case-lambda:
                                              [()               (unsafe-flvector-ref (flstudent-t-sample μ σ ν 1) 0)]
                                              [([n : Integer])  (flvector->list (flstudent-t-sample μ σ ν n))]))
                            (define median (delay ν))
                            (student-t-dist-struct pdf sample cdf inv-cdf -inf.0 +inf.0 median μ σ ν)]))])))


;;;
;;; Tests
;;;

#;(let ()
  (: nearly-equal? : (Real Real Real -> Boolean))
  (define (nearly-equal? eps x y)
    (<= (abs (- x y)) eps))

  (define T (student-t-dist 2))
  (list "Density - PDF"
        (list
         (nearly-equal? (expt 2 -54) (pdf T 1)                                  (/ 1 (* 3 (sqrt 3))))
         (nearly-equal? (expt 2 -54) (cdf T 1)                                  (+ 1/2 (/ 1 (* 2 (sqrt 3)))))
         (nearly-equal? (expt 2 -51) (inv-cdf T (+ 1/2 (/ 1 (* 2 (sqrt 3)))))   1))))

