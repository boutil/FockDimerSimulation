package dimerSim;

import java.io.Serializable;

public class Index implements Serializable{
    final int x;
    final int y;

    public Index(int x, int y) {this.x=x;this.y=y;}

    public boolean isNeighbor (Index other) {
        return (Math.abs(x - other.x) + Math.abs(y - other.y)) == 1;
    }

    public boolean equals (Object o) {
        if (o == this) {
            return true;
        }
        if (!(o instanceof Index)) {
            return false;
        }
        Index c = (Index) o;
        return (c.x == x) & (c.y == y);
    }

    public Index minus(int i, int j) {
        return new Index(x - i, y - j);
    }

    public int l1Dist(int i, int j) {
        return Math.abs(i - x) + Math.abs(j - y);
    }

    public boolean isEven() {
        return ((x + y) % 2) == 0;
    }
}