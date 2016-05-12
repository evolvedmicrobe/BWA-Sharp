using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace Bio.BWA
{
    /// <summary>
    ///  Holds genomic regions based on Binary Search Tree.
    /// </summary>
    public class RegionTree
    {
        #region BTreeNode Class
        [DebuggerDisplay("Value: {Value}")]
        public class BTreeNode
        {
            /// <summary>
            /// Constructor to initialize null node.
            /// </summary>
            public BTreeNode()
            {
                this.Left = this;
                this.Right = this;
            }

            public BTreeNode(Region value, BTreeNode nullNode)
            {
                this.Value = value;
                this.Left = nullNode;
                this.Right = nullNode;
                this.Level = 1;
            }

            public Region Value { get; set; }
            public BTreeNode Left { get; set; }
            public BTreeNode Right { get; set; }
            public int Level { get; set; }

        }
        #endregion

        #region Member variables
        /// <summary>
        /// Holds null node.
        /// </summary>
        private static BTreeNode NullNode = new BTreeNode(new Region(String.Empty, 0, 1), RegionTree.NullNode);

        /// <summary>
        /// Holds Root node.n
        /// </summary>
        private BTreeNode root;

        #endregion

        #region Constructors

        /// <summary>
        /// Initializes an instance of AATree class with specified comparer.
        /// </summary>
        public RegionTree()
        {
            this.root = NullNode;
            this.DefaultValue = default(Region);
        }
        #endregion

        #region Properties
        /// <summary>
        /// Gets or sets Default value.
        /// By default this is set to default value of T.
        /// </summary>
        public Region DefaultValue { get; set; }

        /// <summary>
        /// Gets number of elements present in the AATree.
        /// </summary>
        public long Count { get; private set; }
        #endregion

        #region Methods
        /// <summary>
        /// Tries to add specified value to the AATree.
        /// If the value is already present in the tree then this method returns without adding.
        /// </summary>
        /// <param name="value">Value to add.</param>
        /// <returns>Returns true if value is added successfully, else returns false.</returns>
        public void Add(Region value)
        {
            bool result = false;
            if (this.root == RegionTree.NullNode)
            {
                this.root = new BTreeNode(value, RegionTree.NullNode);
                result = true;
            }

            if (!result)
            {
                BTreeNode node = this.root;
                Stack<BTreeNode> parentNodes = new Stack<BTreeNode>();
                parentNodes.Push (node);
                while (true)
                {
                    int keyCompResult = value.CompareTo(node.Value);
                    if (keyCompResult == 0)
                    {
                        // key already exists.
                        result = false;
                        // If the old node doesn't completely contain the new range
                        // it is possible multiple existing nodes need to be merged, so
                        // we'll look for those. 
                        bool contained = node.Value.FullyContainsRegion (value);
                        node.Value.Merge (value);
                        if (!contained) {
                            // Have to check if we merge any child intervals
                            var additionalNodes = findOverlappingChildNodes (node);
                            foreach (var v in additionalNodes) {
                                node.Value.Merge (v.Item1.Value);
                                DeleteNode (v.Item1, v.Item2);
                            }
                        }
                        break;
                    }
                    else if (keyCompResult < 0)
                    {
                        // go to left.
                        if (node.Left == RegionTree.NullNode)
                        {
                            node.Left = new BTreeNode(value, RegionTree.NullNode);
                            result = true;
                            break;
                        }
                        else
                        {
                            parentNodes.Push(node);
                            node = node.Left;
                        }
                    }
                    else
                    {
                        // go to right.
                        if (node.Right == RegionTree.NullNode)
                        {
                            node.Right = new BTreeNode(value, RegionTree.NullNode);
                            result = true;
                            break;
                        }
                        else
                        {
                            parentNodes.Push(node);
                            node = node.Right;
                        }
                    }
                }

                if (result)
                {
                    parentNodes.Push(node);
                   
                }
            }

            if (result)
            {
                Count++;
            }
        }


        private Tuple<BTreeNode, BTreeNode> FindMinimumAndParent (BTreeNode start, BTreeNode parent)
        {
            if (start.Left == NullNode) {
                return new Tuple<BTreeNode, BTreeNode> (start, parent);
            } else {
                return FindMinimumAndParent (start.Left, start);
            }
        }

        /// <summary>
        /// Deletes the node.
        /// 
        /// http://www.algolist.net/Data_structures/Binary_search_tree/Removal
        /// 
        /// </summary>
        /// <returns>The node.</returns>
        /// <param name="toDelete">To delete.</param>
        /// <param name="parent">Parent.</param>
        private void DeleteNode (BTreeNode toDelete, BTreeNode parent)
        {
            Count--;
            // Is this a leaf node?
            if (toDelete.Left == NullNode && toDelete.Right == NullNode) {
                // just remove it
                if (parent.Left == toDelete) {
                    parent.Left = NullNode;
                } else {
                    parent.Right = NullNode;
                }
            } else if (toDelete.Left == NullNode || toDelete.Right == NullNode) {
                // It has one child
                // Delete the node and replace it with it's child
                var child = toDelete.Left == NullNode ? toDelete.Right : toDelete.Left;
                if (parent.Left == toDelete) {
                    parent.Left = child;
                } else {
                    parent.Right = child;
                }
            } else {
                // It has two children
                // find minimum on the right, we get the parent at the same time
                // to avoid finding it again during a recursion
                var values = FindMinimumAndParent (toDelete.Right, parent);
                var minNode = values.Item1;
                var minsParent = values.Item2;

                // Two cases, the parent is the node we are deleting, or it's something further down
                if (minsParent != toDelete) {
                    // Remove this node
                    Count++; // Account for how we are making a new node somewhere
                    DeleteNode (minNode, minsParent);
                    minNode.Right = toDelete.Right;
                }
                minNode.Left = toDelete.Left;               
                if (parent.Left == toDelete) {
                    parent.Left = minNode;
                } else {
                    parent.Right = minNode;
                }
            }
        }

        private List<Tuple<BTreeNode, BTreeNode>> findOverlappingChildNodes (BTreeNode justMerged) {
            var region = justMerged.Value;
            // List of child/parent pairs
            List<Tuple<BTreeNode, BTreeNode>> overlaps = new List<Tuple<BTreeNode, BTreeNode>> ();
            var toVisit = new Stack<Tuple<BTreeNode, BTreeNode>> ();
            toVisit.Push (new Tuple<BTreeNode, BTreeNode>(justMerged, NullNode));
            while (toVisit.Count != 0) {
                var currentTuple = toVisit.Pop ();
                var currentNode = currentTuple.Item1;
                var currentParent = currentTuple.Item2;
                var comparison = currentNode.Value.CompareTo (region);
                if (comparison == 0 && currentNode != justMerged) {
                    overlaps.Add (currentTuple);
                }
                // No overlap, keep searching in same direction
                if (currentNode.Left != NullNode && comparison >= 0) {
                    toVisit.Push (new Tuple<BTreeNode, BTreeNode>(currentNode.Left, currentNode));
                };
                if (currentNode.Right != NullNode && comparison <= 0) {
                    toVisit.Push (new Tuple<BTreeNode, BTreeNode> (currentNode.Right, currentNode));
                    //toVisit.Push (new Tuple<BTreeNode, BTreeNode>(currentNode.Right, currentNode)); t
                }
            }
            return overlaps;
        }



        /// <summary>
        /// Searches for the specified value in the AATree.
        /// If found returns the node containing the value in node out param, else this param contains NullNode.
        /// </summary>
        /// <param name="value">Value to search.</param>
        /// <param name="node">AATree node.</param>
        /// <returns>Returns true if the value is found, else returns false.</returns>
        internal bool TrySearch(Region value, out Region output)
        {
            bool result = false;
            BTreeNode node;
            node = RegionTree.NullNode;

            var currentNode = this.root;
            while (currentNode != RegionTree.NullNode)
            {
                int compResult = value.CompareTo(currentNode.Value);
                if (compResult == 0)
                {
                    node = currentNode;
                    result = true;
                    break;
                }

                if (compResult < 0)
                {
                    //goto left
                    currentNode = currentNode.Left;
                }
                else
                {
                    //goto right
                    currentNode = currentNode.Right;
                }
            }
            output = node.Value;
            return result;
        }

        /// <summary>
        /// Gets values using InOrder traversal.
        /// </summary>
        public IEnumerable<Region> InOrderTraversal()
        {
            return RegionTree.InOrderTraversal(this.root);
        }

        /// <summary>
        /// Gets values using PreOrder traversal.
        /// </summary>
        public IEnumerable<Region> PreOrderTraversal()
        {
            return RegionTree.PreOrderTraversal(this.root);
        }

        /// <summary>
        /// Gets values using PostOrder traversal.
        /// </summary>
        public IEnumerable<Region> PostOrderTraversal()
        {
            return RegionTree.PostOrderTraversal(this.root);
        }


        // InOrder Traversal implementation.
        private static IEnumerable<Region> InOrderTraversal(BTreeNode node)
        {
            if (node == RegionTree.NullNode)
            {
                yield break;
            }

            Stack<BTreeNode> nodes = new Stack<BTreeNode>();
            BTreeNode currentNode = node;
            while (nodes.Count > 0 || currentNode != RegionTree.NullNode)
            {
                if (currentNode != RegionTree.NullNode)
                {
                    nodes.Push(currentNode);
                    currentNode = currentNode.Left;
                }
                else
                {
                    currentNode = nodes.Pop();
                    yield return currentNode.Value;
                    currentNode = currentNode.Right;
                }
            }
        }

        // PreOrder Traversal implementation.
        private static IEnumerable<Region> PreOrderTraversal(BTreeNode node)
        {
            if (node == RegionTree.NullNode)
            {
                yield break;
            }

            Stack<BTreeNode> stack = new Stack<BTreeNode>();
            stack.Push(node);
            BTreeNode currentNode = RegionTree.NullNode;
            while (stack.Count > 0)
            {
                currentNode = stack.Pop();
                if (currentNode.Right != RegionTree.NullNode)
                {
                    stack.Push(currentNode.Right);
                }

                if (currentNode.Left != RegionTree.NullNode)
                {
                    stack.Push(currentNode.Left);
                }

                yield return currentNode.Value;
            }
        }

        // PostOrder Traversal implementation.
        private static IEnumerable<Region> PostOrderTraversal(BTreeNode node)
        {
            if (node == RegionTree.NullNode)
            {
                yield break;
            }

            Stack<BTreeNode> nodes = new Stack<BTreeNode>();
            BTreeNode currentNode = node;
            while (true)
            {
                if (currentNode != RegionTree.NullNode)
                {
                    if (currentNode.Right != RegionTree.NullNode)
                    {
                        nodes.Push(currentNode.Right);
                    }

                    nodes.Push(currentNode);
                    currentNode = currentNode.Left;
                }
                else
                {
                    if (nodes.Count == 0)
                    {
                        break;
                    }

                    currentNode = nodes.Pop();
                    if (currentNode.Right != RegionTree.NullNode && nodes.Count > 0 && nodes.Peek() == currentNode.Right)
                    {
                        nodes.Pop();  // remove right;
                        nodes.Push(currentNode); // push current
                        currentNode = currentNode.Right;
                    }
                    else
                    {
                        yield return currentNode.Value;
                        currentNode = RegionTree.NullNode;
                    }
                }
            }
        }
        #endregion
    }


}