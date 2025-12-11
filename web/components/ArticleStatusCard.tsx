'use client';

import { useRouter } from 'next/navigation';
import StatusBadge from './StatusBadge';
import { format } from 'date-fns';
import { ja } from 'date-fns/locale/ja';
import { ArticleListItem } from '@/lib/api';

interface ArticleStatusCardProps {
  article: ArticleListItem & {
    pendingApproval?: boolean;
    approvalDeadline?: string;
    kpiSummary?: Record<string, any>;
  };
}

export default function ArticleStatusCard({ article }: ArticleStatusCardProps) {
  const router = useRouter();
  
  const handleClick = () => {
    router.push(`/articles/${article.id}`);
  };
  
  const kpi = article.kpiSummary || {};
  const likesCount = kpi.likes_count || 0;
  const viewsCount = kpi.page_views_count || 0;
  
  return (
    <div 
      onClick={handleClick}
      className={`bg-white dark:bg-gray-800 rounded-lg shadow p-6 hover:shadow-lg transition-shadow cursor-pointer ${
        article.pendingApproval ? 'border-l-4 border-yellow-500' : ''
      }`}
    >
      <div className="flex justify-between items-start mb-4">
        <h3 className="text-lg font-semibold text-gray-900 dark:text-white line-clamp-2">
          {article.title}
        </h3>
        <div className="flex flex-col items-end space-y-1">
          <StatusBadge status={article.uiStatusText} phase={article.phase} />
          {article.pendingApproval && (
            <span className="px-2 py-1 text-xs font-medium bg-yellow-100 text-yellow-800 dark:bg-yellow-900 dark:text-yellow-200 rounded">
              承認待ち
            </span>
          )}
        </div>
      </div>
      <p className="text-sm text-gray-600 dark:text-gray-400 mb-4">
        {article.nextActionHint}
      </p>
      
      {/* パフォーマンスデータ */}
      {(likesCount > 0 || viewsCount > 0) && (
        <div className="flex space-x-4 mb-4 text-xs">
          <div>
            <span className="text-gray-500 dark:text-gray-500">いいね: </span>
            <span className="font-medium text-gray-900 dark:text-white">{likesCount}</span>
          </div>
          <div>
            <span className="text-gray-500 dark:text-gray-500">閲覧: </span>
            <span className="font-medium text-gray-900 dark:text-white">{viewsCount}</span>
          </div>
        </div>
      )}
      
      <div className="flex justify-between items-center text-xs text-gray-500 dark:text-gray-500">
        <span>
          {article.createdAt ? format(new Date(article.createdAt), 'yyyy年MM月dd日 HH:mm', { locale: ja }) : '-'}
        </span>
        <span>
          {article.updatedAt ? format(new Date(article.updatedAt), 'yyyy年MM月dd日 HH:mm', { locale: ja }) : '-'}
        </span>
      </div>
    </div>
  );
}

