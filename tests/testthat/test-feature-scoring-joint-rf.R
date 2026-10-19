.joint_rf_fixture <- function() {
  gp_with_pinned_seed(15L, {
    X <- matrix(rnorm(60 * 10), 60, 10,
                dimnames=list(paste0("L",1:60),paste0("m",1:10)))
    ph <- data.frame(ID=rownames(X), Env=rep(c("E1","E2"),30),
                     T1=3*X[,1]+rnorm(60,sd=.3)+rep(c(0,2),30),
                     T2=2*X[,1]+rnorm(60,sd=.3)+rep(c(0,-3),30))
    list(X=X,ph=ph)
  })
}

test_that("Gaussian matrix targets retain shape without changing single-output targets", {
  Y <- matrix(as.double(1:12),6,2,dimnames=list(NULL,c("T1","T2")))
  expect_identical(gp_py_ml_prepare_target(Y,"gaussian"),Y)
  expect_identical(gp_py_ml_prepare_target(Y[,1,drop=FALSE],"gaussian"),as.numeric(Y[,1]))
  expect_identical(gp_py_ml_prepare_target(1:6,"gaussian"),as.numeric(1:6))
  expect_identical(gp_py_ml_prepare_target(factor(c("b","a")),"multiclass"),c("b","a"))
  expect_error(gp_ml_write_vector_csv(Y,tempfile()),"multi-output|matrix")
})

test_that("importance CLI writes all target columns and retains their record order", {
  d <- .joint_rf_fixture()
  Y <- as.matrix(d$ph[,c("T1","T2")])
  captured <- NULL
  local_mocked_bindings(gp_ml_run_cli=function(command,args,...) {
    y_path <- args[match("--y-train-csv",args)+1L]
    captured <<- as.matrix(read.csv(y_path,check.names=FALSE))
    out <- args[match("--out-dir",args)+1L]
    write.csv(data.frame(importance=seq_len(ncol(d$X))),
              file.path(out,"feature_importance.csv"),row.names=FALSE)
    invisible(NULL)
  },.package="PredictProR")
  imp <- gp_ml_feature_importance_cli("randomforest",d$X,Y,
                                     model_params=list(importance="air",random_state=11L))
  expect_equal(unname(captured),unname(Y),tolerance=1e-12)
  expect_equal(imp,as.numeric(seq_len(ncol(d$X))))
})

test_that("joint RF sends residualized scaled matrix y to AIR and keeps signed scores", {
  d <- .joint_rf_fixture()
  captured <- NULL
  local_mocked_bindings(gp_py_ml_feature_importance=function(model_type,X_train,y_train,
                                                            response_family,model_params,
                                                            class_levels) {
    captured <<- list(X=X_train,y=y_train,params=model_params,family=response_family)
    -seq_len(ncol(X_train))
  },.package="PredictProR")
  meta <- feature_score_predictors(d$X,d$ph,c("T1","T2"),"ID",
                                   scoring_model="RandomForest",joint=TRUE,heter_groups="Env")
  C <- model.matrix(~ Env,d$ph)
  expected <- vapply(c("T1","T2"),function(tr) {
    residual <- lm.fit(C,d$ph[[tr]])$residuals
    residual/sd(residual)
  },numeric(nrow(d$ph)))
  expect_equal(unname(captured$y),unname(expected),tolerance=1e-12)
  expect_equal(dim(captured$X),dim(d$X))
  expect_identical(captured$params$importance,"air")
  expect_identical(captured$family,"gaussian")
  expect_true(all(meta$raw_score<0))
  expect_identical(meta$predictor[1L],"m1")
  expect_identical(unique(meta$scoring_engine),"sklearn RandomForest AIR (multi-output)")
  expect_true(all(is.na(meta$effect_sd)))
  expect_true(all(is.na(meta$chain_nIter)))
})

test_that("RF complete-case drops are warned and environment drops are counted once", {
  d <- .joint_rf_fixture()
  d$ph$T1[1:4] <- NA_real_
  d$ph$Env[4:5] <- NA_character_
  local_mocked_bindings(gp_py_ml_feature_importance=function(X_train,...) {
    seq_len(ncol(X_train))
  },.package="PredictProR")
  expect_warning(meta <- feature_score_predictors(
    d$X,d$ph,c("T1","T2"),"ID",scoring_model="RandomForest",joint=TRUE,heter_groups="Env"),
    "dropped 4 incomplete")
  expect_identical(unique(meta$n_train),55L)
  expect_identical(unique(meta$n_records_dropped),5L)
  expect_identical(unique(meta$n_incomplete_records),4L)
  d <- .joint_rf_fixture()
  d$ph$T1[1:30] <- NA_real_; d$ph$T2[31:60] <- NA_real_
  expect_error(feature_score_predictors(
    d$X,d$ph,c("T1","T2"),"ID",scoring_model="RandomForest",joint=TRUE),"no complete")
})

test_that("joint RF CLI agrees with direct Python AIR on the same matrix targets", {
  py <- skip_if_no_ml_python()
  skip_if_python_lacks(py,c("pandas"),"ML")
  d <- .joint_rf_fixture()
  z <- gp_feature_joint_inputs(d$X,d$ph,c("T1","T2"),"ID","Env",
                               complete_only=TRUE,residualize_markers=FALSE)
  yr <- vapply(seq_len(ncol(z$y)),function(j) lm.fit(z$C,z$y[,j])$residuals,
               numeric(nrow(z$y)))
  params <- list(ntree=120L,mtry=3L,nodesize=2L,n_jobs=1L,
                 random_state=19L,importance="air",air_unranked=0L)
  td <- tempfile("joint_rf_reference_"); dir.create(td)
  on.exit(unlink(td,recursive=TRUE),add=TRUE)
  gp_ml_write_matrix_csv(z$x,file.path(td,"x.csv"))
  gp_ml_write_matrix_csv(yr,file.path(td,"y.csv"))
  gp_ml_write_json(params,file.path(td,"params.json"))
  py_dir <- system.file("python",package="PredictProR")
  script <- c("import sys,json,numpy as np",
              "sys.path.insert(0,sys.argv[1])",
              "from ml_models import _air_importance",
              "from ml_bridge import _read_target_csv",
              "X=np.loadtxt(sys.argv[2],delimiter=',',skiprows=1)",
              "y=_read_target_csv(sys.argv[3],'gaussian')",
              "assert y.ndim==2 and y.shape[1]==2",
              "params=json.load(open(sys.argv[4],encoding='utf-8'))",
              "np.savetxt(sys.argv[5],_air_importance('randomforest','gaussian',X,y,params,None))")
  writeLines(script,file.path(td,"direct.py"))
  status <- system2(py,args=vapply(c(file.path(td,"direct.py"),py_dir,
    file.path(td,"x.csv"),file.path(td,"y.csv"),file.path(td,"params.json"),
    file.path(td,"importance.txt")),shQuote,character(1)),
    stdout=file.path(td,"python.log"),stderr=file.path(td,"python.log"))
  expect_identical(status,0L)
  direct <- scan(file.path(td,"importance.txt"),quiet=TRUE)
  set.seed(50L); before <- .Random.seed
  meta <- feature_score_predictors(
    d$X,d$ph,c("T1","T2"),"ID",joint=TRUE,scoring_model="RandomForest",
    heter_groups="Env",seed=19L,ntree=120L,mtry=3L,nodesize=2L,rf_n_jobs=1L)
  expect_equal(meta$raw_score[match(colnames(d$X),meta$predictor)],direct,tolerance=1e-12)
  expect_identical(meta$predictor[1L],"m1")
  expect_identical(.Random.seed,before)
})
