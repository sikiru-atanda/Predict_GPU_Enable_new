# Regression guard: apply_aliases() deleted identity-mapped settings (e.g.
# "ft_layers"="ft_layers"), so user values for TabTransformer, TabAttention,
# Conv1DNet pooling, NeuralAdditive, MixtureOfExperts and GPNet were silently
# replaced by internal defaults.

test_that("identity aliases keep the user's deep-learning settings", {
  dl <- function(type, args) PredictProR:::gp_dl_bridge_fit_params(type, args)

  ft <- dl("ft_transformer", list(ft_d_model = 64L, ft_heads = 4L, ft_layers = 2L, ft_ff_mult = 2L,
                                  ft_dropout = 0.05, ft_token_dropout = 0, ft_use_cls = FALSE))
  expect_identical(ft$ft_d_model, 64L)
  expect_identical(ft$ft_heads, 4L)
  expect_identical(ft$ft_layers, 2L)
  expect_identical(ft$ft_ff_mult, 2L)
  expect_equal(ft$ft_dropout, 0.05)
  expect_equal(ft$ft_token_dropout, 0)
  expect_false(ft$ft_use_cls)

  saint <- dl("saint", list(saint_d_model = 32L, saint_layers = 1L, saint_heads = 2L))
  expect_identical(saint$saint_d_model, 32L)
  expect_identical(saint$saint_layers, 1L)
  expect_identical(saint$saint_heads, 2L)

  cnn <- dl("cnn", list(cnn_use_max_pool = TRUE, cnn_pool_kernel = 3L, cnn_pool_stride = 3L))
  expect_true(cnn$cnn_use_max_pool)
  expect_identical(cnn$cnn_pool_kernel, 3L)
  expect_identical(cnn$cnn_pool_stride, 3L)

  nam <- dl("nam", list(nam_hidden = c(8L, 4L), nam_l1 = 0.5))
  expect_identical(nam$nam_hidden, c(8L, 4L))
  expect_equal(nam$nam_l1, 0.5)

  moe <- dl("moe", list(moe_n_experts = 3L, moe_temperature = 2, moe_sparse_topk = NA))
  expect_identical(moe$moe_n_experts, 3L)
  expect_equal(moe$moe_temperature, 2)
  # NA = not set: omitted so Python uses its default (was sent as "NA" and crashed)
  expect_false("moe_sparse_topk" %in% names(moe))

  gpdkl <- dl("gp_dkl", list(gp_feature_dim = 16L, gp_num_inducing = 20L))
  expect_identical(gpdkl$gp_feature_dim, 16L)
  expect_identical(gpdkl$gp_num_inducing, 20L)

  # non-identity aliases still map to their canonical names
  mlp <- dl("mlp", list(mlp_neurons_per_layer = c(32L, 16L)))
  expect_identical(mlp$neurons_per_layer, c(32L, 16L))
})
