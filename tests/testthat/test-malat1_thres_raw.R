test_that("malat1_thres_raw returns threshold and cells", {

  data(example_counts)
  data(example_barcodes)

  result <- malat1_thres_raw(
    example_counts,
    example_barcodes,
    print_plots = FALSE
  )

  expect_true("threshold" %in% names(result))
  expect_true("cells" %in% names(result))
  expect_true(length(result$cells) > 0)

})

