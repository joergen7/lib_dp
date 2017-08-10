# lib_dp
###### Dynamic Programming in Erlang.

[![hex.pm](https://img.shields.io/hexpm/v/lib_dp.svg?style=flat-square)](https://hex.pm/packages/lib_dp) [![Build Status](https://travis-ci.org/joergen7/lib_dp.svg?branch=master)](https://travis-ci.org/joergen7/lib_dp)

This library allows the inexact matching of a pair of sequences via Dynamic
Programming. The library supports three different matching schemes:

- global alignment (Needleman-Wunsch)
- global end-space free alignment
- local alignment (Smith-Waterman)

Furthermore, it is possible to provide custom substitution matrices and to set
the score of insertions and deletions. After alignment the similarity score of
both sequences can be queried and different kinds of edit scripts can be
produced.

The algorithm has O(n^2) complexity regarding the length of the input sequences
in both processing time as well as memory consumption.

## Compiling from Source

The following packages should be installed on your machine:

- [git](https://git-scm.com)
- [erlang](http://www.erlang.org/) (OTP 18.0 or higher)
- [rebar3](https://github.com/erlang/rebar3)

Download the git repository and change into the repository's directory:

    git clone https://github.com/joergen7/lib_dp.git
    cd lib_dp

Now, compile the library with rebar3

    rebar3 compile

In addition, you can run unit tests or check the library for discrepancies with
the following two commands:

    rebar3 eunit
    rebar3 dialyzer

You can build the documentation by entering:

    rebar3 edoc

After compilation you can start a shell with the `lib_dp` library available by
entering:

    rebar3 shell

## Example

Start an Erlang shell by entering:

    rebar3 shell

`lib_dp` allows the matching of lists of any given elements. Since Erlang
strings are internally represented as lists of integers, it is straightforward
to match them. Note, however, that the lists must be flat in order for the
matching to be successful.

For this example, we will, however, consider lists of atoms. Let us define two
similar sequences (lists of atoms):

    A = [v, i, n, t, n, e, r].
    B = [w, i, n, t, e, r, s].

First, we initialize a stateful module which holds general parameters about the
way, we want to perform alignments. Here, we will perform a global alignment
(Needleman-Wunsch). Thus we initialize the stateful module by entering:

    Dp = lib_dp:new( global ).

Now, we compile the Dynamic Programming score table from the two strings `A`
and `B`:

    Tbl = Dp:scoretbl( A, B ).

If the score table is small (in this example it is) we can output it on the
console and inspect it:

    lib_dp:print_tbl( Tbl ).

The table will look something like this:

    \ -1 - -2 - -3 - -4 - -5 - -6 - -7 
    | -2 \  0 - -1 - -2 - -3 - -4 - -5 
    | -3 | -1 \  1 -  0 - -1 - -2 - -3 
    | -4 | -2 |  0 \  2 -  1 -  0 - -1 
    | -5 | -3 | -1 |  1 \  1 -  0 - -1 
    | -6 | -4 | -2 |  0 \  2 -  1 -  0 
    | -7 | -5 | -3 | -1 |  1 \  3 -  2 

It displays the score associated with each position pair as well as the
back-links that memorize from which direction the score was generated. Higher
scores mean higher similarity.

For global alignments, the similarity score of both sequences is the value in
the lower right corner of the score table, here: 2. We can extract the
similarity score from a score table by applying the function `get_score`:

    Dp:get_score( Tbl ).

Note that the score will change depending on the alignment strategy.

We can now cunstruct an edit script from the programming table by entering:

    EditScr = Dp:editscr( A, B, Tbl ).

The edit script is a list of pairs in which the first element represents a
symbol from sequence `A` and the second represents a symbol from sequence `B`.
The edit script for this example is the following:

    [{v,w},{i,i},{n,n},{t,t},{n,indel},{e,e},{r,r},{indel,s}]

Eventually, we can simplify the edit script to a summary script, which is just
another way to represent an edit script in memory:

    lib_dp:sumscr( EditScr ).

The resulting summary script looks like this:

    [{mismatch,[v],[w]},
     {match,[i,n,t]},
     {ins,[n]},
     {match,[e,r]},
     {del,[s]}]

## System Requirements

- Erlang OTP 18.0 or higher
- Rebar3 3.0.0 or higher

## Authors

- Jörgen Brandt (joergen7) [joergen.brandt@onlinehome.de](mailto:joergen.brandt@onlinehome.de)

## License

[Apache 2.0](https://www.apache.org/licenses/LICENSE-2.0.html)