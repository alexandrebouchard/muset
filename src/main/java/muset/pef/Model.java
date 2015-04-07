package muset.pef;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import muset.Alphabet;
import muset.Alphabet.Letter;
import muset.Sequence;
import briefj.Indexer;



public final class Model
{
  public final StrTaxonSuffStat stSuffStat;
  public final Alphabet enc;
  public final int nStates;
  public final Indexer<Object> stateIndexer;
  public final int startState, endState;
  public final int epsilon;
  public final int boundaryId;
  public final Letter BOUNDARY_SYMBOL; // used to mark the end of the strings
  public Model(StrTaxonSuffStat stSuffStat, Alphabet enc, Indexer<Object> stateIndexer, int startState, int endState)
  {
    this.stSuffStat = stSuffStat;
    this.enc = ensureContainsBoundarySymbol(enc);
    this.stateIndexer = stateIndexer;
    this.nStates = stateIndexer.size();
    this.startState = startState;
    this.endState = endState;
    this.epsilon = this.enc.indexer.size();
    this.BOUNDARY_SYMBOL = enc.getLetter(_BOUNDARY_SYMBOL_STR);
    this.boundaryId = enc.indexer.o2i(BOUNDARY_SYMBOL);
  }
  
  private static final String _BOUNDARY_SYMBOL_STR = "#";
  private static Alphabet ensureContainsBoundarySymbol(
      Alphabet alphabet)
  {
    alphabet = new Alphabet(alphabet);
    if (!alphabet.containsLetter(_BOUNDARY_SYMBOL_STR))
      alphabet.getLetter(_BOUNDARY_SYMBOL_STR);//indexer.addToIndex(BOUNDARY_SYMBOL);
    return alphabet;
  }
  public int epsilon() { return epsilon; }
  public List<Input> allInputs()
  {
    List<Input> result = new ArrayList<Input>();
    for (int s1 = 0; s1 < nStates; s1++)
      for (int stss = 0; stss < stSuffStat.valuesIndexer.size(); stss++)
        result.add(new Input(s1, stss, this));
    return result;
  }
  public List<Output> allOutputs()
  {
    List<Output> result = new ArrayList<Output>();
    for (int s1 = 0; s1 < nStates; s1++)
      for (int sym1 = 0; sym1< enc.indexer.size() + 1; sym1++)
        for (int sym2 = 0; sym2 < enc.indexer.size() + 1; sym2++)
          if (sym1 != epsilon() || sym2 != epsilon())
            result.add(new Output(s1, sym1, sym2, this));
    return result;
  }
  
  public static final String SUB = "SUB";
  public static final String INS = "INS";
  public static final String DEL = "DEL";
  public static Model stdModel(Alphabet alphabet) { return stdBranchSpecificModel(alphabet, null); }
  public static Model stdBranchSpecificModel(Alphabet alphabet, StrTaxonSuffStat stss)
  {
    Indexer<Object> stateIndexer = new Indexer<Object>();
    stateIndexer.addToIndex(SUB);
    stateIndexer.addToIndex(INS);
    stateIndexer.addToIndex(DEL);
    int specialState = stateIndexer.o2i(SUB);
    return new Model(stss,
        alphabet, stateIndexer, specialState, specialState);
  }
  
  public int charIdAt(Sequence str, int xpos, int dx)
  {
    if (dx == 0) return epsilon();
    else if (dx == 1) 
    {
      int result = enc.indexer.o2i(str.letterAt(xpos));
      if (result == -1)
        throw new RuntimeException("Unknown character: " + str.letterAt(xpos) +  " in sequence: " + str);
      return result;
    }
    else throw new RuntimeException();
  }
  
  public static class ThreeStatesBaseMeasure implements bayonet.regression.BaseMeasures<Input,Output>
  {
    private static final long serialVersionUID = 1L;
    //    private final Model model;
    private final SortedSet<Output> 
     from01 = new TreeSet<Output>(), 
     from2 = new TreeSet<Output>();
    private final int s0,s1,s2;
    public ThreeStatesBaseMeasure(Model model)
    {
      int bound = model.boundaryId;
      s0 = model.stateIndexer.o2i(SUB);
      s1 = model.stateIndexer.o2i(INS);
      s2 = model.stateIndexer.o2i(DEL);
      if (model.stateIndexer.size() != 3) throw new RuntimeException();
//      this.model = model;
      // from 0/1:
      {
        for (int top = 0; top < model.enc.indexer.size(); top++)
          for (int bot = 0; bot < model.enc.indexer.size(); bot++)
            if (bothOrNeitherBound(bound, top, bot))
              from01.add(new Output(s0, top, bot,model));
        for (int bot = 0; bot < model.enc.indexer.size(); bot++)
          if (bot != bound)
            from01.add(new Output(s1, model.epsilon(), bot, model));
        for (int top = 0; top < model.enc.indexer.size(); top++)
          if (top !=bound)
            from01.add(new Output(s2, top, model.epsilon(), model));
      }
      // from 2:
      {
        for (int top = 0; top < model.enc.indexer.size(); top++)
          for (int bot = 0; bot < model.enc.indexer.size(); bot++)
            if (bothOrNeitherBound(bound, top, bot))
              from2.add(new Output(s0, top, bot,model));
        for (int top = 0; top < model.enc.indexer.size(); top++)
          if (top != bound)
            from2.add(new Output(s2, top, model.epsilon(), model));
      }
    }
    private boolean bothOrNeitherBound(int bound, int top, int bot)
    {
      if (top == bound && bot == bound) return true;
      return top != bound && bot != bound;
    }
    @Override
    public SortedSet<Output> support(Input input)
    {
      if (input.state1 == s0 || input.state1 == s1)
        return from01;
      else if (input.state1 == s2)
        return from2;
      else
        throw new RuntimeException();
    }
  }
  
}
