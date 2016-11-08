%% -*- erlang -*-
%%
%% Dynamic Programming in Erlang
%%
%% Copyright 2016 Jörgen Brandt
%%
%% Licensed under the Apache License, Version 2.0 (the "License");
%% you may not use this file except in compliance with the License.
%% You may obtain a copy of the License at
%%
%%    http://www.apache.org/licenses/LICENSE-2.0
%%
%% Unless required by applicable law or agreed to in writing, software
%% distributed under the License is distributed on an "AS IS" BASIS,
%% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%% See the License for the specific language governing permissions and
%% limitations under the License.

%% @author Jörgen Brandt <brandjoe@hu-berlin.de>

%% @doc Dynamic Programming library.

-module( lib_dp ).
-author( "Jorgen Brandt <brandjoe@hu-berlin.de>" ).

-export( [new/1, new/3, scoretbl/3, editscr/4, sumscr/1] ).

-export( [find_max_score/3, extend_backward/4, extend_forward/4,
          find_local_matches/4] ).

-ifdef( TEST ).
-include_lib( "eunit/include/eunit.hrl" ).
-endif.


%% ============================================================================
%% Record definitions
%% ============================================================================

-record( lib_dp, {strategy, subst_fun, indel_score} ).

%% ============================================================================
%% Types
%% ============================================================================

% The score table is a hash map from a pair of positions in two sequences to
% a pair holding a score and a back-link.
-type score_table()    :: #{{pos_integer(), pos_integer()} => score_lnk()}.

% A sequence of deletions as appears in a summary script.
-type del()            :: {del, [_]}.

% An edit script is a list of pairs. The elements of a pair are either the atom
% `indel' or a token from the corresponding sequence.
-type edit_script()    :: [{_, _}].

% A sequence of insertions as appears in a summary script.
-type ins()            :: {ins, [_]}.

% A sequence of matches as appears in a summary script.
-type match()          :: {match, [_]}.

% A sequence of mismatches as appears in a summary script.
-type mismatch()       :: {mismatch, [_], [_]}.

% A pair of positions in two sequences
-type pos_pair()       :: {non_neg_integer(), non_neg_integer()}.

% A score link pair consists of a score and a back-link which is either a
% position pair or the atom `undef'.
-type score_lnk()      :: {number(), pos_pair() | undef}.

% A strategy with which to build the Dynamic Programming score table. Global
% alignment (Needleman-Wunsch), global end-space free alignment, or local
% alignment (Smith-Waterman) are possible options.
-type strategy()       :: global | global_endfree | local.

% A substitution function which given a pair of tokens returns the score that is
% associated with substituting the first for the latter.
-type subst_fun()      :: fun( ( _, _ ) -> number() ).

% A condensed edit script.
-type summary_script() :: [match() | mismatch() | ins() | del()].


%% ============================================================================
%% API functions
%% ============================================================================


%% @doc Constructs an edit script from two sequences and the Dynamic Programming
%%      score table.

-spec editscr( TLst1, TLst2, Tbl, ModArg ) -> edit_script()
when TLst1  :: [_],
     TLst2  :: [_],
     Tbl    :: score_table(),
     ModArg :: #lib_dp{}.

editscr( TLst1, TLst2, Tbl, ModArg=#lib_dp{ strategy=global } ) ->
  ILen = length( TLst1 ),
  JLen = length( TLst2 ),
  RLst1 = lists:reverse( TLst1 ),
  RLst2 = lists:reverse( TLst2 ),
  {Scr, {[], []}} = editscr_fixpnt( {0, 0}, {ILen, JLen}, {RLst1, RLst2}, Tbl,
                                    [], ModArg ),
  Scr;

editscr( TLst1, TLst2, Tbl, ModArg=#lib_dp{ strategy=global_endfree } ) ->
  ILen = length( TLst1 ),
  JLen = length( TLst2 ),
  {{VIndex, JLen}, VScore} = find_vertical_end( {ILen, JLen}, Tbl ),
  {{ILen, HIndex}, HScore} = find_horizontal_end( {ILen, JLen}, Tbl ),
  if
    VScore > HScore ->
      Suffix = lists:nthtail( VIndex, TLst1 ),
      RPrefix = lists:reverse( lists:sublist( TLst1, VIndex ) ),
      RLst2 = lists:reverse( TLst2 ),
      {Scr, {[], []}} = editscr_fixpnt( {0, 0}, {VIndex, JLen},
                                        {RPrefix, RLst2}, Tbl, [], ModArg ),
      Scr++[{T, indel} || T <- Suffix];
    true ->
      Suffix = lists:nthtail( HIndex, TLst2 ),
      RPrefix = lists:reverse( lists:sublist( TLst2, HIndex ) ),
      RLst1 = lists:reverse( TLst1 ),
      {Scr, {[], []}} = editscr_fixpnt( {0, 0}, {ILen, HIndex},
                                        {RLst1, RPrefix}, Tbl, [], ModArg ),
      Scr++[{indel, T} || T <- Suffix]
  end;

editscr( S1, S2, Tbl, ModArg=#lib_dp{ strategy=local } ) ->
  S1Len = length( S1 ),
  S2Len = length( S2 ),
  RMatches = lists:reverse( find_local_matches( {0, 0}, {S1Len, S2Len}, Tbl,
                                                ModArg ) ),
  RS1 = lists:reverse( S1 ),
  RS2 = lists:reverse( S2 ),

  F = fun

        % If we are at position {0, 0} just return the accumulator.
        % Note that we can be sure that both reverse sequences are empty
        % and the reverse list of local matches is also empty.
        F( {0, 0}, {[], []}, [], Acc ) ->
          Acc;

        % If the reverse list of local matches is empty but we haven't reached
        % the position {0, 0} yet then we move to that position.
        F( PosPair, RSeqPair, [], Acc )    ->
          {Scr, RSeqPair1} = editscr_fixpnt( {0, 0}, PosPair, RSeqPair, Tbl,
                                             [], ModArg ),
          F( {0, 0}, RSeqPair1, [], Scr++Acc );

        % If we are at a position, equal to the end position of the next match
        % then we process that match.        
        F( PosPair, RSeqPair, [{MatchStart, PosPair}|T], Acc ) ->
          {Scr, RSeqPair1} = editscr_fixpnt( MatchStart, PosPair, RSeqPair,
                                             Tbl, [], ModArg ),
          F( MatchStart, RSeqPair1, T, Scr++Acc );

        % If we are at a position, different from the end position of the next
        % match then we move to that end position.
        F( PosPair, RSeqPair, Matches=[{_, MatchEnd}|_], Acc ) ->
          {Scr, RSeqPair1} = editscr_fixpnt( MatchEnd, PosPair, RSeqPair, Tbl,
                                             [], ModArg ),
          F( MatchEnd, RSeqPair1, Matches, Scr++Acc )
      end,

  F( {S1Len, S2Len}, {RS1, RS2}, RMatches, [] ).


%% @doc A simplified version of the function `new/3' to construct a stateful
%%      module holding the parameters for the Dynamic Programming scheme.
%%
%%      Only a strategy needs to be provided being either the symbol `global',
%%      `global_endfree', or `local'.

-spec new( Strategy :: strategy() ) -> #lib_dp{}.

new( Strategy ) ->

  SubstFun = fun( A, A ) -> 1;
                ( _, _ ) -> -1
             end,

  new( Strategy, SubstFun, -1 ).


%% @doc Constructs a new stateful module holding the parameters for the Dynamic
%%      Programming scheme.
%%
%%      The stateful module is initialized with
%%      <ul>
%%        <li>
%%          A strategy, being either `global' for global alignment
%%          (Needleman-Wunsch), `global_endfree' for global end-space free
%%          alignment, or `local' for local alignment (Smith-Waterman).
%%        </li>
%%        <li>
%%          A substitution function which, given two symbols returns the score
%%          of exchanging them. Typically, a positive score is assigned to
%%          matching symbols and a negative score is assigned to non-matching
%%          symbols.
%%        </li>
%%        <li>
%%          A score for making an insertion or deletion. This should also be a
%%          negative score.
%%        </li>
%%      </ul>
%%
%%      A simpler initialization function where only the the strategy has to be
%%      provided is {@link new/1}.
%%
%%      The returned stateful module is used with the functions
%%      {@link scoretbl/3} to create a score table from two sequences and
%%      {@link editscr/4} to  construct an edit script from a score table.

-spec new( Strategy, SubstFun, Indel ) -> #lib_dp{}
when Strategy :: strategy(),
     SubstFun :: subst_fun(),
     Indel    :: number().

new( Strategy, SubstFun, Indel ) ->
  {?MODULE, Strategy, SubstFun, Indel}.


%% @doc Creates the Dynamic Programming score table from two sequences.

scoretbl( TLst1, TLst2, ModArg ) ->

  F = fun( T1, {Tbl1, I} ) ->

        G = fun( T2, {Tbl2, J} ) ->
              D = calc_scorelnk( {I, J}, {T1, T2}, Tbl2, ModArg ),
              {Tbl2#{ {I, J} => D }, J+1}
            end,

        {Tbl3, _} = lists:foldl( G, {Tbl1, 1}, TLst2 ),
        {Tbl3, I+1}
      end,

  {Tbl4, _} = lists:foldl( F, {#{}, 1}, TLst1 ),
  Tbl4.


%% @doc Creates a summary script from an edit script.

-spec sumscr( edit_script() ) -> summary_script().

sumscr( [] ) ->
  [];

sumscr( [{indel, Y}|T] ) ->
  sumscr_fixpnt( T, {del, [Y]}, [] );

sumscr( [{X, indel}|T] ) ->
  sumscr_fixpnt( T, {ins, [X]}, [] );

sumscr( [{X, X}|T] ) ->
  sumscr_fixpnt( T, {match, [X]}, [] );

sumscr( [{X, Y}|T] ) ->
  sumscr_fixpnt( T, {mismatch, [X], [Y]}, [] ).


%% ============================================================================
%% Internal functions
%% ============================================================================

%% @doc The fixpoint function that constructs the edit script from the Dynamic
%%      Programming score table.
%%
%%      The function constructs a list of pairs, the left element being a
%%      symbol from the first sequence or the atom `indel' and the right
%%      element being a symbol from the second sequence or the atom `indel'.
%%
%%      The function takes six arguments:
%%      <ul>
%%        <li>
%%          A pair of positions in both sequences at which the edit script
%%          starts. The edit script exludes the actual positions. Hence, to
%%          have an edit script start at the beginning of both sequences, the
%%          start position pair would be `{0, 0}'.
%%        </li>
%%        <li>
%%          A pair of positions in both sequences at which the edit script
%%          ends. The edit script includes the actual positions. Hence, to
%%          have an edit script end at the end of two sequences S1 and S2 the
%%          end position pair would be `{length( S1 ), length( S2 )}'.
%%        </li>
%%        <li>
%%          A pair of reversed sequences. The sequences need to be in reversed
%%          form because the construction of the edit script starts from the
%%          end of both sequences. Hence for two sequences S1 and S2 the
%%          reverse sequence pair would be
%%          `{lists:reverse( S1 ), lists:reverse( S2 )}'.
%%        </li>
%%        <li>
%%          The Dynamic Programming score table constructed for the two
%%          sequences.
%%        </li>
%%        <li>
%%          The edit script accumulator which should be the empty list `[]'.
%%        </li>
%%        <li>
%%          The stateful module holding the strategy with which the Dynamic
%%          Programming score table was constructed.
%%        </li>
%%      </ul>

-spec editscr_fixpnt( StartPair, EndPair, RSeqPair, Tbl, Acc, ModArg ) ->
  {edit_script(), {[_], [_]}}
when StartPair :: {non_neg_integer(), non_neg_integer()},
     EndPair   :: {non_neg_integer(), non_neg_integer()},
     RSeqPair  :: {[_], [_]},
     Tbl       :: score_table(),
     Acc       :: edit_script(),
     ModArg    :: #lib_dp{}.

editscr_fixpnt( {I0, J0}, {I0, J0}, RSeqPair, _, Acc, _ ) ->
  {Acc, RSeqPair};

editscr_fixpnt( Start={I0, J0}, {I0, J}, {L1, [H2|T2]}, Tbl, Acc, ModArg ) ->
  editscr_fixpnt( Start, {I0, J-1}, {L1, T2}, Tbl, [{indel, H2}|Acc], ModArg );


editscr_fixpnt( Start={I0, J0}, {I, J0}, {[H1|T1], L2}, Tbl, Acc, ModArg ) ->
  editscr_fixpnt( Start, {I-1, J0}, {T1, L2}, Tbl, [{H1, indel}|Acc], ModArg );

editscr_fixpnt( Start={I0, _}, Pos={I, J}, {L1=[H1|T1], L2=[H2|T2]}, Tbl, Acc,
                ModArg ) ->

  {_, Lnk} = get_scorelnk( Pos, Tbl, ModArg ),

  case Lnk of
    undef ->
      if
        I > I0 ->
          editscr_fixpnt( Start, {I-1, J}, {T1, L2}, Tbl, [{H1, indel}|Acc],
                          ModArg );
        true ->
          editscr_fixpnt( Start, {I, J-1}, {L1, T2}, Tbl, [{indel, H2}|Acc],
                          ModArg )
      end;
    {ILnk, JLnk} ->
      if
        I > ILnk ->
          if
            J > JLnk ->
              editscr_fixpnt( Start, {I-1, J-1}, {T1, T2}, Tbl, [{H1, H2}|Acc],
                              ModArg );
            true  ->
              editscr_fixpnt( Start, {I-1, J}, {T1, L2}, Tbl, [{H1, indel}|Acc],
                              ModArg )
          end;
        true ->
          editscr_fixpnt( Start, {I, J-1}, {L1, T2}, Tbl, [{indel, H2}|Acc],
                          ModArg )
      end
  end.

%% @doc Extends a local alignment from a given position in backward direction.

extend_backward( Pos={I, _}, {I, _}, _, _ ) -> Pos;
extend_backward( Pos={_, J}, {_, J}, _, _ ) -> Pos;

extend_backward( Pos, ScopeStart, Tbl, ModArg ) ->

  {_, Lnk} = get_scorelnk( Pos, Tbl, ModArg ),

  if
    Lnk =:= undef -> Pos;
    true          -> extend_backward( Lnk, ScopeStart, Tbl, ModArg )
  end.

%% @doc Extends a local alignment from a given position in forward direction.

extend_forward( Pos={_, JEnd}, {_, JEnd}, _, _ ) -> Pos;
extend_forward( Pos={IEnd, _}, {IEnd, _}, _, _ ) -> Pos;

extend_forward( Pos={I, J}, ScopeEnd, Tbl, ModArg ) ->

  {A, _} = get_scorelnk( {I+1, J}, Tbl, ModArg ),
  {B, _} = get_scorelnk( {I, J+1}, Tbl, ModArg ),
  {C, _} = get_scorelnk( {I+1, J+1}, Tbl, ModArg ),

  case lists:max( [A, B, C] ) > 0 of
    false -> Pos;  
    true  ->
      if
        C > A ->
          if
            C > B ->
              extend_forward( {I+1, J+1}, ScopeEnd, Tbl, ModArg );
            true  ->
              extend_forward( {I, J+1}, ScopeEnd, Tbl, ModArg )
          end;
        true ->
          if
            A > B ->
              extend_forward( {I+1, J}, ScopeEnd, Tbl, ModArg );
            true  ->
              extend_forward( {I, J+1}, ScopeEnd, Tbl, ModArg )
          end
      end
  end.


%% @doc Finds the optimal end in the second sequence in a global end-space free
%%      alignment.
%%
%%      The function produces a pair consisting of a position pair and its
%%      score.

-spec find_horizontal_end( EndPair, Tbl ) -> {pos_pair(), number()}
when EndPair :: pos_pair(),
     Tbl     :: score_table().

find_horizontal_end( {ILen, JLen}, Tbl ) ->
  find_max_score( {ILen-1, 0}, {ILen, JLen}, Tbl ).


find_local_matches( {I, _}, {I, _}, _, _ ) -> [];
find_local_matches( {_, J}, {_, J}, _, _ ) -> [];

find_local_matches( ScopeStart, ScopeEnd, Tbl, ModArg ) ->


  {HitCenter, Score} = find_max_score( ScopeStart, ScopeEnd, Tbl ),

  if
    Score =< 0 -> [];
    true       ->


      HitStart = extend_backward( HitCenter, ScopeStart, Tbl, ModArg ),
      HitEnd = extend_forward( HitCenter, ScopeEnd, Tbl, ModArg ),
      
      PreMatchLst = find_local_matches( ScopeStart, HitStart, Tbl, ModArg ),
      PostMatchLst = find_local_matches( HitEnd, ScopeEnd, Tbl, ModArg ),

      PreMatchLst++[{HitStart, HitEnd}|PostMatchLst]
  end.


%% @doc Finds the position in a Dynamic Programming score table with the highest
%%      score.

-spec find_max_score( ScopeStart, ScopeEnd, Tbl ) -> {pos_pair(), number()}
when ScopeStart :: pos_pair(),
     ScopeEnd   :: pos_pair(),
     Tbl        :: score_table().

find_max_score( {I0, J0}, {I1, J1}, Tbl ) ->

  FilterPred = fun( {I, J}, _ ) 
               when I > I0 andalso
                    J > J0 andalso
                    I =< I1 andalso
                    J =< J1 -> true;
                  ( _, _ )  -> false
               end,


  Acc = fun( Pos, {Score, _}, {_, ScoreMax} ) when Score > ScoreMax ->
             {Pos, Score};
           ( _, _, Acc ) ->
             Acc
        end,

  Tbl1 = maps:filter( FilterPred, Tbl ),

  [Pos|_] = maps:keys( Tbl1 ),
  {Score, _} = maps:get( Pos, Tbl1 ),

  maps:fold( Acc, {Pos, Score}, Tbl1 ).





%% @doc Finds the optimal end in the first sequence in a global end-space free
%%      alignment.
%%
%%      The function produces a pair consisting of a position pair and its
%%      score.

-spec find_vertical_end( EndPair, Tbl ) -> {pos_pair(), number()}
when EndPair :: pos_pair(),
     Tbl     :: score_table().

find_vertical_end( {ILen, JLen}, Tbl ) ->
  find_max_score( {0, JLen-1}, {ILen, JLen}, Tbl ).
  

%% @doc Looks up the score of an alignment at a given position.

-spec get_scorelnk( Pos, Tbl, ModArg ) -> score_lnk()
when Pos    :: pos_pair(),
     Tbl    :: score_table(),
     ModArg :: #lib_dp{}.

get_scorelnk( {0, 0}, _, _ ) ->
  {0, undef};

get_scorelnk( {I, 0}, _, #lib_dp{ strategy=global, indel_score=IndelScore } ) ->
 {I*IndelScore, {I-1, 0}};

get_scorelnk( {0, J}, _, #lib_dp{ strategy=global, indel_score=IndelScore } ) ->
  {J*IndelScore, {0, J-1}};

get_scorelnk( {I, 0}, _, _ ) ->
  {0, {I-1, 0}};

get_scorelnk( {0, J}, _, _ ) ->
  {0, {0, J-1}};

get_scorelnk( Key, Tbl, _ ) ->
  maps:get( Key, Tbl ).

%% @doc Calculates the score of an alignment at a given position.
%%
%%      The function produces a number representing the score at the given
%%      position depending on the score of preceding positions and the given
%%      alignment strategy.
%%
%%      The four arguments are:
%%      <ul>
%%        <li>
%%          A pair of positions in both sequences.
%%        </li>
%%        <li>
%%          A pair of tokens corresponding to the position pair.
%%        </li>
%%        <li>
%%          The currently available Dynamic Programming score table that
%%          memoizes previously calculated scores.
%%        </li>
%%        <li>
%%          A stateful module holding the alignment strategy, substitution
%%          function, and indel score.
%%        </li>
%%      </ul>
%%
%%      Note that this score calculation function does not work with an
%%      insufficiently filled score table.
%%
%%      The function is used to create the Dynamic Programming score table in
%%      {@link scoretbl/3}.

-spec calc_scorelnk( Pos, TokenPair, Tbl, ModArg ) -> score_lnk()
when Pos       :: pos_pair(),
     TokenPair :: {_, _},
     Tbl       :: score_table(),
     ModArg    :: #lib_dp{}.

calc_scorelnk( {I, J},
               {T1, T2},
               Tbl,
               ModArg=#lib_dp{ strategy    = Strategy,
                               subst_fun   = SubstFun,
                               indel_score = IndelScore } ) ->

  {A0, _} = get_scorelnk( {I-1, J}, Tbl, ModArg ),
  {B0, _} = get_scorelnk( {I, J-1}, Tbl, ModArg ),
  {C0, _} = get_scorelnk( {I-1, J-1}, Tbl, ModArg ),

  A1 = A0+IndelScore,
  B1 = B0+IndelScore,
  C1 = C0+SubstFun( T1, T2 ),

  case Strategy =:= local andalso lists:max( [A1, B1, C1] ) =< 0 of
    true  -> {0, undef};
    false ->
      if
        C1 > A1 ->
          if
            C1 > B1 -> {C1, {I-1, J-1}};
            true    -> {B1, {I, J-1}}
          end;
        true    ->
          if
            A1 > B1 -> {A1, {I-1, J}};
            true    -> {B1, {I, J-1}}
          end
      end
  end.



% termination
sumscr_fixpnt( [], {Op, Acc}, MajAcc ) ->
  lists:reverse( [{Op, lists:reverse( Acc )}|MajAcc] );

sumscr_fixpnt( [], {mismatch, A, B}, MajAcc ) ->
  lists:reverse( [{mismatch, lists:reverse( A ), lists:reverse( B )}|MajAcc] );

% from match to X
sumscr_fixpnt( [{indel, Y}|T], {match, Acc}, MajAcc ) ->
  sumscr_fixpnt( T, {del, [Y]}, [{match, lists:reverse( Acc )}|MajAcc] );

sumscr_fixpnt( [{X, indel}|T], {match, Acc}, MajAcc ) ->
  sumscr_fixpnt( T, {ins, [X]}, [{match, lists:reverse( Acc )}|MajAcc] );

sumscr_fixpnt( [{X, X}|T], {match, Acc}, MajAcc ) ->
  sumscr_fixpnt( T, {match, [X|Acc]}, MajAcc );

sumscr_fixpnt( [{X, Y}|T], {match, Acc}, MajAcc ) ->
  sumscr_fixpnt( T, {mismatch, [X], [Y]}, [{match, lists:reverse( Acc )}|MajAcc] );

% from insert to X
sumscr_fixpnt( [{indel, Y}|T], {ins, Acc}, MajAcc ) ->
  sumscr_fixpnt( T, {del, [Y]}, [{ins, lists:reverse( Acc )}|MajAcc] );

sumscr_fixpnt( [{X, indel}|T], {ins, Acc}, MajAcc ) ->
  sumscr_fixpnt( T, {ins, [X|Acc]}, MajAcc );

sumscr_fixpnt( [{X, X}|T], {ins, Acc}, MajAcc ) ->
  sumscr_fixpnt( T, {match, [X]}, [{ins, lists:reverse( Acc )}|MajAcc] );

sumscr_fixpnt( [{X, Y}|T], {ins, Acc}, MajAcc ) ->
  sumscr_fixpnt( T, {mismatch, [X], [Y]}, [{ins, lists:reverse( Acc )}|MajAcc] );

% from mismatch to X
sumscr_fixpnt( [{indel, Y}|T], {mismatch, A, B}, MajAcc ) ->
  sumscr_fixpnt( T, {del, [Y]}, [{mismatch, lists:reverse( A ), lists:reverse( B )}|MajAcc] );

sumscr_fixpnt( [{X, indel}|T], {mismatch, A, B}, MajAcc ) ->
  sumscr_fixpnt( T, {ins, [X]}, [{mismatch, lists:reverse( A ), lists:reverse( B )}|MajAcc] );

sumscr_fixpnt( [{X, X}|T], {mismatch, A, B}, MajAcc ) ->
  sumscr_fixpnt( T, {match, [X]}, [{mismatch, lists:reverse( A ), lists:reverse( B )}|MajAcc] );

sumscr_fixpnt( [{X, Y}|T], {mismatch, A, B}, MajAcc ) ->
  sumscr_fixpnt( T, {mismatch, [X|A], [Y|B]}, MajAcc );

% from del to X
sumscr_fixpnt( [{indel, Y}|T], {del, Acc}, MajAcc ) ->
  sumscr_fixpnt( T, {del, [Y|Acc]}, MajAcc );

sumscr_fixpnt( [{X, indel}|T], {del, Acc}, MajAcc ) ->
  sumscr_fixpnt( T, {ins, [X]}, [{del, lists:reverse( Acc )}|MajAcc] );

sumscr_fixpnt( [{X, X}|T], {del, Acc}, MajAcc ) ->
  sumscr_fixpnt( T, {match, [X]}, [{del, lists:reverse( Acc )}|MajAcc] );

sumscr_fixpnt( [{X, Y}|T], {del, Acc}, MajAcc ) ->
  sumscr_fixpnt( T, {mismatch, [X], [Y]}, [{del, lists:reverse( Acc )}|MajAcc] ).




%% ============================================================================
%% Unit tests
%% ============================================================================

-ifdef( TEST ).

find_vertical_end_finds_max1_test() ->
  Tbl = #{
    {1, 15} => -3,
    {2, 15} => -4,
    {3, 15} => 20,
    {4, 15} => -1
  },
  ?assertEqual( {{3, 15}, 20}, find_vertical_end( {4, 15}, Tbl ) ).

find_vertical_end_finds_max2_test() ->
  Tbl = #{
    {1, 15} => -3,
    {2, 15} => -4,
    {3, 15} => -1,
    {4, 15} => 20
  },
  ?assertEqual( {{4, 15}, 20}, find_vertical_end( {4, 15}, Tbl ) ).

find_vertical_end_finds_max3_test() ->
  Tbl = #{
    {1, 15} => 20,
    {2, 15} => -4,
    {3, 15} => -3,
    {4, 15} => -1
  },
  ?assertEqual( {{1, 15}, 20}, find_vertical_end( {4, 15}, Tbl ) ).


find_vertical_end_returns_zero_on_empty_s1_test() ->
  ?assertEqual( {{0, 0}, 0}, find_vertical_end( {4, 15}, #{} ) ).

find_max_score_finds_max_score_test() ->
  Tbl = #{
    {1, 1} => -3,
    {2, 1} => -4,
    {3, 1} => -3,
    {4, 1} => -1,
    {1, 2} => -5,
    {2, 2} => 20,
    {3, 2} => -3,
    {4, 2} => -1
  },
  ?assertEqual( {{2, 2}, 20}, find_max_score( {0, 0}, {4, 2}, Tbl ) ).
-endif.


