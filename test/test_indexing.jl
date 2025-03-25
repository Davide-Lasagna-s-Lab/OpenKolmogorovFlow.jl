@testset "_reindex                               " begin
    # test data m = 1
    #  0 0 0
    #  1 0 0
    #  2 0 0
    # -1 0 0
    @test OpenKolmogorovFlow._reindex( 0, 0, 1) == (0, 0)
    @test OpenKolmogorovFlow._reindex( 0, 1, 1) == (0, 1)
    @test OpenKolmogorovFlow._reindex( 1, 1, 1) == (1, 1)
    @test OpenKolmogorovFlow._reindex( 1, 0, 1) == (1, 0)
    @test OpenKolmogorovFlow._reindex(-1, 1, 1) == (3, 1)
    @test OpenKolmogorovFlow._reindex(-1, 0, 1) == (3, 0)

    # test data m = 2
    #  0 0 0 0
    #  1 0 0 0
    #  2 0 0 0
    #  3 0 0 0
    # -2 0 0 0
    # -1 0 0 0
    @test OpenKolmogorovFlow._reindex( 0, 0, 2) == (0, 0)
    @test OpenKolmogorovFlow._reindex( 1, 0, 2) == (1, 0)
    @test OpenKolmogorovFlow._reindex( 2, 0, 2) == (2, 0)
    @test OpenKolmogorovFlow._reindex(-2, 0, 2) == (4, 0)
    @test OpenKolmogorovFlow._reindex(-1, 0, 2) == (5, 0)
end

# @testset "js and ks                              " begin
#     @test OpenKolmogorovFlow.js(4) == (0, 1, 2)
#     @test OpenKolmogorovFlow.js(6) == (0, 1, 2,  3)
#     @test OpenKolmogorovFlow.ks(4) == (0, 1, 2, -1)
#     @test OpenKolmogorovFlow.ks(6) == (0, 1, 2, 3, -2, -1)
# end

# @testset "JtoI                                   " begin
#     @test OpenKolmogorovFlow.JtoI( 3) == 4
#     @test OpenKolmogorovFlow.JtoI( 0) == 1
#     @test OpenKolmogorovFlow.JtoI(-3) == 4
#     @test OpenKolmogorovFlow.KtoI( 0, 4)  == 1
#     @test OpenKolmogorovFlow.KtoI( 1, 4)  == 2
#     @test OpenKolmogorovFlow.KtoI(-1, 4)  == 4
#     @test OpenKolmogorovFlow.KtoI(-1, 6)  == 6
#     @test OpenKolmogorovFlow.KtoI(-2, 6)  == 5
# end

# @testset "ItoK                                   " begin
#     # these are the wave numbers for k
#     data_k = [0  0  0
#               1  1  1
#               2  2  2
#              -1 -1 -1]
#     for i in 1:12
#         @test data_k[i] == OpenKolmogorovFlow.ItoK(i, 4)
#     end
# end

# @testset "ItoJ                                   " begin
#     # these are the wave numbers for j
#     data_j = [0  1  2
#               0  1  2
#               0  1  2
#               0  1  2]
#     for i in 1:12
#         @test data_j[i] == OpenKolmogorovFlow.ItoJ(i, 4)
#     end
# end

# @testset "KJtoI                                  " begin
#     # data = [1+2im  9+10im 0+0im
#     #         3+4im 11+12im 0+0im
#     #         0+0im  0+0im  0+0im
#     #         3-4im 15+16im 0+0im]
#     @test OpenKolmogorovFlow.KJtoI( 0,  0, 1) ==  1
#     @test OpenKolmogorovFlow.KJtoI( 1,  0, 1) ==  2
#     @test OpenKolmogorovFlow.KJtoI(-1,  0, 1) ==  4
    
#     @test OpenKolmogorovFlow.KJtoI( 0, -1, 1) ==  5
#     @test OpenKolmogorovFlow.KJtoI( 0,  1, 1) ==  5
    
#     @test OpenKolmogorovFlow.KJtoI( 1,  1, 1) ==  6
#     @test OpenKolmogorovFlow.KJtoI(-1, -1, 1) ==  6

#     @test OpenKolmogorovFlow.KJtoI( 1, -1, 1) ==  8
#     @test OpenKolmogorovFlow.KJtoI(-1,  1, 1) ==  8
# end

# @testset "rectify                                " begin
#     @test OpenKolmogorovFlow.rectify(1,   1)     == 1
#     @test OpenKolmogorovFlow.rectify(1,   0)     == 1
#     @test OpenKolmogorovFlow.rectify(1,  -1)     == 1
#     @test OpenKolmogorovFlow.rectify(1 + im,  1) == 1 + im
#     @test OpenKolmogorovFlow.rectify(1 + im,  0) == 1 + im
#     @test OpenKolmogorovFlow.rectify(1 + im, -1) == 1 - im
# end