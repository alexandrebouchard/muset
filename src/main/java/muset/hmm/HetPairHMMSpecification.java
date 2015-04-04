package muset.hmm;

/**
 * Convention: start state is 0, final state is nStates() - 1
 * @author bouchard
 *
 */
public interface HetPairHMMSpecification
{
  /**
   * x and deltaX are positions in the top string (str1)
   * y and deltaY are positions in the bot string (str2)
   * For a substitution, deltaX=deltaY=1, and the character being substituted is [x,x+deltaX)=str1.charAt(x), same idea for y
   * For an insert (i.e. insertion in str2), deltaX=0, deltaY=1, the character being inserted in [y,y+deltaY)=str2.charAt(y)
   * Same idea for delete, i.e. insert in str1
   * Assumes for now that for all nonzero weight transitions, at most one of deltaX, deltaY is nonzero, i.e. that the 
   * pure epsilon transitions have been removed
   */
  public double logWeight(int prevState, int currentState, int x, int y, int deltaX, int deltaY);
  public int nStates();
  public int startState();
  public int endState();
}