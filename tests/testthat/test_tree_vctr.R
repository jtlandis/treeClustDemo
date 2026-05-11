set.seed(123)

x <- matrix(rnorm(20), 5, 4)
y <- matrix(rnorm(12, 4), 3, 4)
z <- rbind(x, y)

xtree <- tree_vctr(hclust(dist(x)))
ytree <- tree_vctr(hclust(dist(y)))
merged_tree <- c(xtree, ytree)

test_that("individual trees have distinct factor levels", {
  expect_identical(
    tree_id(xtree),
    rep(1L, 5)
  )
  expect_identical(
    tree_id(ytree),
    rep(1L, 3)
  )
  expect_false(
    levels(vctrs::field(xtree, "which_tree")) ==
      levels(vctrs::field(ytree, "which_tree"))
  )
})

test_that("merged contains", {
  expect_identical(
    tree_id(merged_tree),
    c(rep(1L, 5), rep(2L, 3))
  )
  expect_identical(
    node_id(encode_obsv(merged_tree)),
    c(seq_len(5), seq_len(3))
  )
  smpl <- sample(merged_tree)
  expect_identical(
    levels(vctrs::field(smpl, "which_tree")),
    levels(vctrs::field(merged_tree, "which_tree"))
  )
  expect_identical(
    trees(merged_tree),
    trees(smpl)
  )
})

test_that("subsetting removes unneeded trees", {
  tree_0 <- merged_tree[1:3]
  expect_identical(
    tree_id(tree_0),
    rep(1L, 3)
  )
  expect_true(
    length(trees(tree_0)) == 1
  )
  expect_true(
    length(levels(vctrs::field(tree_0, "which_tree"))) == 1L
  )
})

test_that("subset can match original", {
  tree_0 <- merged_tree[1:3]
  expect_identical(
    vctrs::vec_match(tree_0, merged_tree),
    1:3
  )
  tree_1 <- merged_tree[c(6, 8)]
  expect_identical(
    vctrs::vec_match(tree_1, merged_tree),
    c(6L, 8L)
  )
})


test_that("ptypes", {
  # expect_identical(
  #   vec_ptype(merged_tree),
  #   vec_slice(merged_tree, 0L)
  # )

  expect_identical(
    vctrs::vec_ptype2(xtree, ytree),
    vctrs::vec_ptype(merged_tree)
  )
})

test_that("descendants are different based on encoding", {
  ord <- xtree
  obs <- encode_obsv(xtree)
  desc_ord <- descendants(inner(ord))
  desc_obs <- descendants(inner(obs))
  expect_false(identical(desc_ord, desc_obs))
  # however if matched against the orignal node id
  expect_true(
    identical(
      lapply(desc_ord, match, table = node_id(ord)),
      lapply(desc_obs, match, table = node_id(obs))
    )
  )
})

test_that("generate_inner_slice yields the same slices regardless of encoding", {
  expect_identical(
    generate_inner_slice(xtree)$children,
    generate_inner_slice(encode_obsv(xtree))$children
  )
  expect_identical(
    generate_inner_slice(merged_tree)$children,
    generate_inner_slice(encode_obsv(merged_tree))$children
  )
})
