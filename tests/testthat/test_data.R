context('data')

test_that("test='HM27.address' gives correct data", {
    sesameDataCache("HM27.address")
    dt <- sesameDataGet('HM27.address')
    expect_is(dt, "list")
})

